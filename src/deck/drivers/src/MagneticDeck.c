#define DEBUG_MODULE "MagneticDeck"
#include "debug.h"
#include "deck.h"
#include "deck_analog.h"
#include "stm32f4xx.h"
#include "stm32f4xx_gpio.h"
#include "stm32f4xx_adc.h"
#include "stm32f4xx_dma.h"

uint32_t DMA_Buffer[19];

static void magneticInit()
{
    ADC_TypeDef *ADC_n = ADC1;
    DMA_Stream_TypeDef *DMA_Stream = DMA2_Stream4;
    uint32_t DMA_Channel = DMA_Channel_0;
    IRQn_Type DMA_stream = DMA2_Stream4_IRQn;
    // RCC_APB2_Peripherals RCC_APB2_Peripheral = RCC_APB2Periph_ADC1;
    // gpio init
    GPIO_init(GPIO_PinSource3);
    // DMA init
    DMA_inititalization(DMA_Stream, DMA_Buffer, ADC_n, DMA_Channel, DMA_stream, 19);
    // adc init
    ADC_init_DMA_mode(RCC_APB2Periph_ADC1, ADC_n);
    // dma start
    // adc start

    // Call the ADC_DMA_start function
    ADC_DMA_start(ADC_n, ADC_Channel_7, 1, ADC_SampleTime_15Cycles);

    DEBUG_PRINT("MAGNETIC init ended!\n");
}

static bool magneticTest()
{
    DEBUG_PRINT("MAGNETIC test passed!\n");
    return true;
}

static const DeckDriver magneticDriver = {
    .name = "Magnetic",

    // .usedGpio = DECK_USING_PA3, // controlla questo

    .init = magneticInit,
    .test = magneticTest,
};

DECK_DRIVER(magneticDriver);