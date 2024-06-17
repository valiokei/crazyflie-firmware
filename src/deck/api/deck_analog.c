/*
 *    ||          ____  _ __
 * +------+      / __ )(_) /_______________ _____  ___
 * | 0xBC |     / __  / / __/ ___/ ___/ __ `/_  / / _ \
 * +------+    / /_/ / / /_/ /__/ /  / /_/ / / /_/  __/
 *  ||  ||    /_____/_/\__/\___/_/   \__,_/ /___/\___/
 *
 * Crazyflie control firmware
 *
 * Copyright (C) 2015 Bitcraze AB
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, in version 3.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 *
 * deck_analog.c - Arduino-compatible analog input implementation
 */

#include "deck.h"
#include "debug.h"
#include "stm32fxxx.h"
#include "stm32f4xx_dma.h"
#include "stm32f4xx_misc.h"
#include "stm32f4xx_adc.h"
#include "stm32f4xx_rcc.h"
#include "stm32f4xx_gpio.h"
#include "arm_math.h"

static uint32_t stregResolution;
static uint32_t adcRange;

void adcInit(void)
{
  /*
   * Note: This function initializes only ADC2, and only for single channel, single conversion mode. No DMA, no interrupts, no bells or whistles.
   */

  /* Note that this de-initializes registers for all ADCs (ADCx) */
  ADC_DeInit();

  /* Define ADC init structures */
  ADC_InitTypeDef ADC_InitStructure;
  ADC_CommonInitTypeDef ADC_CommonInitStructure;

  /* Populates structures with reset values */
  ADC_StructInit(&ADC_InitStructure);
  ADC_CommonStructInit(&ADC_CommonInitStructure);

  /* enable ADC clock */
  RCC_APB2PeriphClockCmd(RCC_APB2Periph_ADC2, ENABLE);

  /* init ADCs in independent mode, div clock by two */
  ADC_CommonInitStructure.ADC_Mode = ADC_Mode_Independent;
  ADC_CommonInitStructure.ADC_Prescaler = ADC_Prescaler_Div6; /* HCLK = 168MHz, PCLK2 = 84MHz, ADCCLK = 42MHz (when using ADC_Prescaler_Div2) */
  ADC_CommonInitStructure.ADC_DMAAccessMode = ADC_DMAAccessMode_Disabled;
  ADC_CommonInitStructure.ADC_TwoSamplingDelay = ADC_TwoSamplingDelay_15Cycles;
  // ADC_CommonInitStructure.ADC_TwoSamplingDelay = ADC_TwoSamplingDelay_5Cycles;
  ADC_CommonInit(&ADC_CommonInitStructure);

  /* Init ADC2: 12bit, single-conversion. For Arduino compatibility set 10bit */
  analogReadResolution(12);

  /* Enable ADC2 */
  ADC_Cmd(ADC2, ENABLE);
}

static uint16_t analogReadChannel(uint8_t channel)
{
  /* According to datasheet, minimum sampling time for 12-bit conversion is 15 cycles. */
  ADC_RegularChannelConfig(ADC2, channel, 1, ADC_SampleTime_15Cycles);

  /* Start the conversion */
  ADC_SoftwareStartConv(ADC2);

  /* Wait until conversion completion */
  while (ADC_GetFlagStatus(ADC2, ADC_FLAG_EOC) == RESET)
    ;

  /* Get the conversion value */
  return ADC_GetConversionValue(ADC2);
}

uint16_t analogRead(const deckPin_t pin)
{
  assert_param(deckGPIOMapping[pin.id].adcCh > -1);

  /* Now set the GPIO pin to analog mode. */

  /* Enable clock for the peripheral of the pin.*/
  RCC_AHB1PeriphClockCmd(deckGPIOMapping[pin.id].periph, ENABLE);

  /* Populate structure with RESET values. */
  GPIO_InitTypeDef GPIO_InitStructure;
  GPIO_StructInit(&GPIO_InitStructure);

  /* Initialise the GPIO pin to analog mode. */
  GPIO_InitStructure.GPIO_Pin = deckGPIOMapping[pin.id].pin;
  GPIO_InitStructure.GPIO_Mode = GPIO_Mode_AN;
  GPIO_InitStructure.GPIO_Speed = GPIO_Speed_25MHz;
  GPIO_InitStructure.GPIO_PuPd = GPIO_PuPd_NOPULL;

  /* TODO: Any settling time before we can do ADC after init on the GPIO pin? */
  GPIO_Init(deckGPIOMapping[pin.id].port, &GPIO_InitStructure);

  /* Read the appropriate ADC channel. */
  return analogReadChannel((uint8_t)deckGPIOMapping[pin.id].adcCh);
}

void analogReference(uint8_t type)
{
  /*
   * TODO: We should probably support the Arduino EXTERNAL type here.
   * TODO: Figure out which voltage reference to compensate with.
   */
  assert_param(type == 0 /* DEFAULT */);
}

void analogReadResolution(uint8_t bits)
{
  ADC_InitTypeDef ADC_InitStructure;

  assert_param((bits >= 6) && (bits <= 12));

  adcRange = 1 << bits;
  switch (bits)
  {
  case 12:
    stregResolution = ADC_Resolution_12b;
    break;
  case 10:
    stregResolution = ADC_Resolution_10b;
    break;
  case 8:
    stregResolution = ADC_Resolution_8b;
    break;
  case 6:
    stregResolution = ADC_Resolution_6b;
    break;
  default:
    stregResolution = ADC_Resolution_12b;
    break;
  }

  /* Init ADC2 witch new resolution */
  ADC_InitStructure.ADC_Resolution = stregResolution;
  ADC_InitStructure.ADC_ScanConvMode = DISABLE;
  ADC_InitStructure.ADC_ContinuousConvMode = DISABLE;
  ADC_InitStructure.ADC_ExternalTrigConvEdge = 0;
  ADC_InitStructure.ADC_ExternalTrigConv = 0;
  ADC_InitStructure.ADC_DataAlign = ADC_DataAlign_Right;
  ADC_InitStructure.ADC_NbrOfConversion = 1;
  ADC_Init(ADC2, &ADC_InitStructure);
}

float analogReadVoltage(const deckPin_t pin)
{
  float voltage;

  voltage = analogRead(pin) * VREF / adcRange;

  return voltage;
}

// DMA ADC SECTION
// inspired by https://community.st.com/t5/stm32-mcus-products/stm32f4-discovery-adc-dma-double-buffer/td-p/422343
// and https://forum.bitcraze.io/viewtopic.php?t=2598&start=10

// #define ADC1_DR ((uint32_t)0x4001244C) #arraySize 20000 __IO uint16_t DMA_Buffer[arraySize];
// #define DMA_Str DMA2_Stream4
void GPIO_init_analog(const deckPin_t pin)
{
  assert_param(deckGPIOMapping[pin.id].adcCh > -1);
  RCC_AHB1PeriphClockCmd(deckGPIOMapping[pin.id].periph, ENABLE);

  /* Populate structure with RESET values. */
  GPIO_InitTypeDef GPIO_InitStructure;
  GPIO_StructInit(&GPIO_InitStructure);
  /* Initialise the GPIO pin to analog mode. */
  DEBUG_PRINT("pin.id: %d\n", pin.id);
  GPIO_InitStructure.GPIO_Pin = deckGPIOMapping[pin.id].pin;
  GPIO_InitStructure.GPIO_Mode = GPIO_Mode_AN;
  // GPIO_InitStructure.GPIO_Speed = GPIO_Speed_25MHz; // non so se è necessario
  GPIO_InitStructure.GPIO_PuPd = GPIO_PuPd_NOPULL;
  GPIO_Init(deckGPIOMapping[pin.id].port, &GPIO_InitStructure);
}

void setADC_Off(ADC_TypeDef *ADCx)
{
  // https: // www.st.com/resource/en/reference_manual/rm0090-stm32f405415-stm32f407417-stm32f427437-and-stm32f429439-advanced-armbased-32bit-mcus-stmicroelectronics.pdf
  //   ADC on-off control
  // The ADC is powered on by setting the ADON bit in the ADC_CR2 register. When the ADON
  // bit is set for the first time, it wakes up the ADC from the Power-down mode.
  // Conversion starts when either the SWSTART or the JSWSTART bit is set.
  // You can stop conversion and put the ADC in power down mode by clearing the ADON bit. In
  // this mode the ADC consumes almost no power (only a few µA).
  return;
}

void ADC_init_DMA_mode(uint32_t RCC_APB2Periph_ADCx, ADC_TypeDef *ADCx)
{
  // Penso che questo debba essere fatto dopo l'init del dma

  // Check the ADC number to be initialized
  // should be ADC1

  /* Note that this de-initializes registers for all ADCs (ADCx) */
  ADC_DeInit();

  /* Define ADC init structures */
  ADC_InitTypeDef ADC_InitStructure;
  ADC_CommonInitTypeDef ADC_CommonInitStructure;

  /* Populates structures with reset values */
  ADC_StructInit(&ADC_InitStructure);
  ADC_CommonStructInit(&ADC_CommonInitStructure);

  /* enable ADC clock */
  RCC_APB2PeriphClockCmd(RCC_APB2Periph_ADCx, ENABLE);

  /* enable ADC in indipendent mode for a initial init*/
  // TODO: check the prescaler for the speed
  ADC_CommonInitStructure.ADC_Mode = ADC_Mode_Independent;
  ADC_CommonInitStructure.ADC_Prescaler = ADC_Prescaler_Div6; /* HCLK = 168MHz, PCLK2 = 84MHz, ADCCLK = 42MHz (when using ADC_Prescaler_Div2) */
  ADC_CommonInitStructure.ADC_DMAAccessMode = ADC_DMAAccessMode_Disabled;
  // ADC_CommonInitStructure.ADC_TwoSamplingDelay = ADC_TwoSamplingDelay_5Cycles;
  ADC_CommonInit(&ADC_CommonInitStructure);

  /* init DMA mode ADC setting*/
  // TODO: check the paramenters values
  // parametri ispirati a https://community.st.com/t5/stm32-mcus-products/stm32f4-discovery-adc-dma-double-buffer/td-p/422343#:~:text=I%20post%20them,should%20be%20illustrative
  ADC_InitStructure.ADC_ContinuousConvMode = ENABLE;
  ADC_InitStructure.ADC_Resolution = ADC_Resolution_12b;
  ADC_InitStructure.ADC_DataAlign = ADC_DataAlign_Right;
  ADC_InitStructure.ADC_ExternalTrigConvEdge = ADC_ExternalTrigConvEdge_None; // TODO: controlla questo parametro

  // ADC_InitStructure.ADC_ExternalTrigConvEdge = ADC_ExternalTrigConvEdge_Rising; // TODO: controlla questo parametro
  ADC_ITConfig(ADCx, ADC_IT_OVR, ENABLE);
  ADC_InitStructure.ADC_NbrOfConversion = 1;
  ADC_InitStructure.ADC_ScanConvMode = DISABLE;
  ADC_Init(ADCx, &ADC_InitStructure);
}

void ADC_DMA_start(ADC_TypeDef *ADC_n, uint8_t ADC_Channel, uint8_t Rank, uint8_t ADC_SampleTime)
{
  // DEBUG_PRINT("****************** ADC_DMA_start initialized !****************** \n");
  /* According to datasheet, minimum sampling time for 12-bit conversion is 15 cycles. */
  // TODO: check the correct sampling time to insert
  //    questo preso da sito stm esempio ispirazione
  // ADC_RegularChannelConfig(ADC1, ADC_Channel_7, 1, ADC_SampleTime_480Cycles);
  //    questo preso da funzione analogReadChannel()
  // ADC_RegularChannelConfig(ADC2, channel, 1, ADC_SampleTime_15Cycles);
  ADC_ContinuousModeCmd(ADC_n, ENABLE);
  ADC_RegularChannelConfig(ADC_n, ADC_Channel, Rank, ADC_SampleTime);

  /*enabling the DMA mode for regular channels group*/
  ADC_DMACmd(ADC_n, ENABLE);

  /*enabling the generation of DMA requests continuously at the end
           of the last DMA transfer*/
  ADC_DMARequestAfterLastTransferCmd(ADC_n, ENABLE);

  // Enable ADC DMA
  ADC_Cmd(ADC_n, ENABLE);
  // NVIC_EnableIRQ(ADC_IRQn);

  // test = 30;
  // Start ADC Conversion
  ADC_SoftwareStartConv(ADC_n);
  // DEBUG_PRINT("****************** ADC_DMA_start ended !****************** \n");
}

void DMA_IRQ_enable(DMA_Stream_TypeDef *DMA_Stream, IRQn_Type DMA_IRQ)
{
  // DMA_ITConfig(DMA_Stream, DMA_IT_HT | DMA_IT_TC, ENABLE);
  DMA_ITConfig(DMA_Stream, DMA_IT_TC | DMA_IT_TE, ENABLE);

  // Enable DMA1 channel IRQ Channel
  NVIC_InitTypeDef NVIC_InitStructure;
  NVIC_InitStructure.NVIC_IRQChannel = DMA_IRQ;
  NVIC_InitStructure.NVIC_IRQChannelPreemptionPriority = 14; // questi sono i valori che indicano la priorita
  NVIC_InitStructure.NVIC_IRQChannelSubPriority = 14;
  NVIC_InitStructure.NVIC_IRQChannelCmd = ENABLE;
  NVIC_Init(&NVIC_InitStructure);
}

void DMA_inititalization(uint32_t RCC_DMA_Peripheral, DMA_Stream_TypeDef *DMA_Stream, uint32_t *DMA_Buffer, ADC_TypeDef *ADC_n, uint32_t DMA_Channel, IRQn_Type DMA_IRQ, uint16_t BufferSize)
{
  RCC_AHB1PeriphClockCmd(RCC_DMA_Peripheral, ENABLE);
  DMA_InitTypeDef DMA_InitStructure;
  DMA_StructInit(&DMA_InitStructure);

  //==Configure DMA2 - Stream 4
  DMA_DeInit(DMA_Stream); // Set DMA registers to default values
  DMA_InitStructure.DMA_Channel = DMA_Channel;
  DMA_InitStructure.DMA_PeripheralBaseAddr = (uint32_t)&ADC_n->DR;  // Source address
  DMA_InitStructure.DMA_Memory0BaseAddr = (uint32_t)&DMA_Buffer[0]; // Destination address
  DMA_InitStructure.DMA_BufferSize = BufferSize;                    // Set the buffer size

  DMA_InitStructure.DMA_PeripheralInc = DMA_PeripheralInc_Disable;
  DMA_InitStructure.DMA_MemoryInc = DMA_MemoryInc_Enable;
  DMA_InitStructure.DMA_PeripheralDataSize = DMA_PeripheralDataSize_Word; // source size - 32bit
  DMA_InitStructure.DMA_MemoryDataSize = DMA_MemoryDataSize_Word;         // destination size = 32b
  DMA_InitStructure.DMA_Mode = DMA_Mode_Normal;
  DMA_InitStructure.DMA_Priority = DMA_Priority_Low;
  DMA_InitStructure.DMA_FIFOMode = DMA_FIFOMode_Disable;
  DMA_InitStructure.DMA_FIFOThreshold = DMA_FIFOThreshold_HalfFull;
  DMA_InitStructure.DMA_MemoryBurst = DMA_MemoryBurst_Single;
  DMA_InitStructure.DMA_PeripheralBurst = DMA_PeripheralBurst_Single;
  DMA_Init(DMA_Stream, &DMA_InitStructure); // Initialize the DMA

  DMA_Cmd(DMA_Stream, ENABLE); // Enable the DMA - Stream x
  DMA_IRQ_enable(DMA_Stream, DMA_IRQ);
}
