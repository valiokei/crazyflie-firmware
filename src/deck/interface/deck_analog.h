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
 * deck_analog.h - Arduino-compatible analog input API
 */

#ifndef __DECK_ANALOG_H__
#define __DECK_ANALOG_H__

#include <stdint.h>
#include "deck_constants.h"
#include "arm_math.h"

/* Voltage reference types for the analogReference() function. */
#define DEFAULT 0
#define VREF 3.0

void adcInit(void);

uint16_t analogRead(const deckPin_t pin);

void analogReference(uint8_t type);

void analogReadResolution(uint8_t bits);

/*
 * Read the voltage on a deck pin.
 * @param[in] pin   deck pin to measure.
 * @return          voltage in volts
 */
float analogReadVoltage(const deckPin_t pin);

// DMA ADC SECTION
void GPIO_init_analog(const deckPin_t pin);

void setADC_Off(ADC_TypeDef *ADCx);

void ADC_init_DMA_mode(uint32_t RCC_APB2Periph_ADCx, ADC_TypeDef *ADCx);

void ADC_DMA_start(ADC_TypeDef *ADC_n, uint8_t ADC_Channel, uint8_t Rank, uint8_t ADC_SampleTime);

void DMA_IRQ_enable(DMA_Stream_TypeDef *DMA_Stream, IRQn_Type DMA_IRQ);

void DMA_inititalization(uint32_t RCC_DMA_Peripheral, DMA_Stream_TypeDef *DMA_Stream, uint32_t *DMA_Buffer, ADC_TypeDef *ADC_n, uint32_t DMA_Channel, IRQn_Type DMA_IRQ, uint16_t BufferSize);

#endif
