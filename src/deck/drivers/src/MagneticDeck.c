#define DEBUG_MODULE "MagneticDeck"
#include "debug.h"
#include "deck.h"
#include "deck_analog.h"
#include "stm32f4xx.h"
#include "stm32f4xx_gpio.h"
#include "stm32f4xx_adc.h"
#include "stm32f4xx_dma.h"
#include "stm32f4xx_rcc.h"
#include "FreeRTOS.h"
#include "semphr.h"
#include "stm32f4xx_misc.h"
#include "timers.h"
#include "stm32f4xx_tim.h"
#include "arm_math.h"
#include "FlattopWinFromPython.h"
#include "arm_const_structs.h"
#include "task.h"
#include "log.h"
#include "param.h"

#define CONFIG_DEBUG = y

// ADC DMA configuration
#define ARRAY_SIZE 2048
uint16_t DMA_Buffer[ARRAY_SIZE];
// float32_t adc_buf_float[ARRAY_SIZE];
ADC_TypeDef *ADC_n = ADC1;
DMA_Stream_TypeDef *DMA_Stream = DMA2_Stream4;
uint32_t DMA_Channel = DMA_Channel_0;
IRQn_Type DMA_IRQ = DMA2_Stream4_IRQn;

// #define DMA_IRQ DMA2_Stream0_IRQn
uint8_t ADC_Channel = ADC_Channel_3;
volatile uint8_t test = 0;
static bool isInit = false;

volatile uint8_t ADC_Done = 0;

volatile uint8_t ADC_Status = 0;
volatile uint8_t DMA_status = 0;
volatile uint8_t DMA_IRQ_status = 0;

volatile uint16_t firstValue = 0;
volatile uint16_t FirstVolt = 0;
volatile float32_t logNeroAmpl = 0;

uint32_t maxindex;
float32_t maxval;
uint32_t valueToExclude = 5;

// FFT parameters
#define FFT_SIZE ARRAY_SIZE
float32_t fft_input[FFT_SIZE];
float32_t fft_output[FFT_SIZE];
float32_t fft_magnitude[FFT_SIZE / 2];
arm_rfft_fast_instance_f32 fft_instance;
uint32_t fft_length = FFT_SIZE;

// 2^12 DOVE 12 SONO I BIT DELL'ADC DEL MCU = 4096
#define ADC_LEVELS 4096
#define ADC_MAX_VOLTAGE 3.3

// Resonance Freqs Anchors in Hz
#define NeroResFreq 117000
// #define NeroIdx 321
#define NeroIdx 2
#define GialloResFreq 97000
#define GialloIdx 1325
#define GrigioResFreq 107000
#define GrigioIdx 1462
#define RossoResFreq 87000
#define RossoIdx 1189

static uint32_t IDx = NeroIdx;

// #define SAMPLES 128
// const int samplingFreq = 1000;
// float32_t mySine[SAMPLES];
// float32_t rfft_output[SAMPLES];
// float32_t test_output[SAMPLES / 2];
// float32_t maxvalue = 10000;
// uint32_t maxindex = 100;
// void performFakeFFT()
// {
//     int peakFrequency;
//     int f = 150;
//     int A = 58;
//     // Create input signal
//     for (int i = 0; i < SAMPLES; i++)
//     {
//         mySine[i] = A * arm_sin_f32(2 * PI * f * i / samplingFreq);
//     }
//     // Create RFFT instance
//     arm_rfft_fast_instance_f32 S;
//     arm_rfft_fast_init_f32(&S, SAMPLES);

//     // RFFT transform
//     arm_rfft_fast_f32(&S, mySine, rfft_output, 0);
//     // Calculate magnitude of imaginary coefficients
//     arm_cmplx_mag_f32(rfft_output, test_output, SAMPLES / 2);
//     // Set DC component to 0
//     test_output[0] = 0;
//     // Obtain peak frequency
//     arm_max_f32(test_output, SAMPLES / 2, &maxvalue, &maxindex);

//     peakFrequency = maxindex * samplingFreq / SAMPLES;

//     DEBUG_PRINT("Peak frequency %d \n\r", peakFrequency);

//     DEBUG_PRINT("Max Value:[%ld]:%f Output=[", maxindex, 2 * maxvalue / SAMPLES);
//     DEBUG_PRINT("]\r\n");
// }

// funzione che fa l'fft e prende in input il puntatore ad un buffer di float
void performFFT(float32_t *Input_buffer_pointer, float32_t *Output_buffer_pointer)
{
    // DEBUG_PRINT("performFFT started!\n");
    // perform FFT on first half of buffer
    arm_q15_to_float((q15_t *)Input_buffer_pointer, Output_buffer_pointer, FFT_SIZE);
    int p = 0;
    for (p = 0; p < FFT_SIZE; p++)
    {
        // Fattore di scala della conversione della funzione, guarda doc per capire
        Output_buffer_pointer[p] = Output_buffer_pointer[p] * 32768;

        // ADCLevelsToVolt
        Output_buffer_pointer[p] = Output_buffer_pointer[p] * (ADC_MAX_VOLTAGE / ADC_LEVELS);
    }
    // applico la finestratura flattop
    arm_mult_f32(Output_buffer_pointer, (float32_t *)flattop_2048_lut, Output_buffer_pointer, FFT_SIZE);

    // eseguo FFT
    arm_rfft_fast_f32(&fft_instance, Output_buffer_pointer, fft_output, 0);

    // normalizzo l'output
    // forse c'e una funzione per fare questo in automatico
    int i = 0;
    for (i = 0; i < FFT_SIZE; i++)
    {
        fft_output[i] = fft_output[i] / FFT_SIZE;
    }

    // calcolo le ampiezze della fft
    arm_cmplx_mag_f32(fft_output, fft_magnitude, FFT_SIZE / 2);

    // leggo le ampiezze per ciascun ancora

    // computing the index of the resonance frequency

    arm_max_f32(&fft_magnitude[valueToExclude], (FFT_SIZE / 2) - valueToExclude, &maxval, &maxindex);

    IDx = maxindex + valueToExclude;
    float32_t NeroAmpl = fft_magnitude[IDx];
    // DEBUG_PRINT("NeroAmpl: %f\n", NeroAmpl);

    // float32_t GialloAmpl = fft_magnitude[GialloIdx];
    // float32_t GrigioAmpl = fft_magnitude[GrigioIdx];
    // float32_t RossoAmpl = fft_magnitude[RossoIdx];

    // compute the log10 of the amplitude
    float32_t logNeroAmpl = pow(10, log10(NeroAmpl));
    DEBUG_PRINT("NeroAmpl_log10: %f\n", logNeroAmpl);

    // TODO CAPIRE QUI COME MANDARE I DATI

    // if ((options->combinedAnchorPositionOk || options->anchorPosition[current_anchor].timestamp) &&
    //     (diff < (OUTLIER_TH * stddev)))
    // {
    //     distanceMeasurement_t dist;
    //     dist.distance = state.distance[current_anchor];
    //     dist.x = options->anchorPosition[current_anchor].x;
    //     dist.y = options->anchorPosition[current_anchor].y;
    //     dist.z = options->anchorPosition[current_anchor].z;
    //     dist.anchorId = current_anchor;
    //     dist.stdDev = 0.25;
    //     estimatorEnqueueDistance(&dist);
    // }
}

static void mytask(void *param)
{
    DEBUG_PRINT("Wait for system starting\n");
    systemWaitStart();
    DEBUG_PRINT("System Started\n");

    // gpio init
    GPIO_init(DECK_GPIO_RX2);
    // DMA init
    DMA_inititalization(RCC_AHB1Periph_DMA2, DMA_Stream, DMA_Buffer, ADC_n, DMA_Channel, DMA_IRQ, ARRAY_SIZE);
    // adc init
    ADC_init_DMA_mode(RCC_APB2Periph_ADC1, ADC_n);
    // Call the ADC_DMA_start function

    ADC_DMA_start(ADC_n, ADC_Channel, 1, ADC_SampleTime_15Cycles);

    arm_rfft_fast_init_f32(&fft_instance, fft_length);
    // DEBUG_PRINT("FFT initialized\n");
    while (1)
    {
        // DEBUG_PRINT("While\n");
        if (ADC_Done == 1)
        {

            performFFT(DMA_Buffer, fft_input);
            firstValue = DMA_Buffer[1000];
            FirstVolt = firstValue * ADC_MAX_VOLTAGE / ADC_LEVELS;
            // DEBUG_PRINT("First Value: %d\n", firstValue);

            DMA_inititalization(RCC_AHB1Periph_DMA2, DMA_Stream, DMA_Buffer, ADC_n, DMA_Channel, DMA_IRQ, ARRAY_SIZE);
            ADC_init_DMA_mode(RCC_APB2Periph_ADC1, ADC_n);
            ADC_DMA_start(ADC_n, ADC_Channel, 1, ADC_SampleTime_15Cycles);

            ADC_Done = 0;
        }
        else
        {
            // DEBUG_PRINT("ADC_Done is 0\n");
        }

        vTaskDelay(M2T(10));
    }
}

void DMA2_Stream4_IRQHandler(void)
{
    // DEBUG_PRINT("DMA2_Stream4_IRQHandler\n");
    // if (DMA_GetITStatus(DMA_Stream, DMA_IT_HTIF0)) //&& (xSemaphoreTake(semaphoreHalfBuffer, 0) == pdTRUE))
    // {
    //     // ADC_Cmd(ADC_n, DISABLE);
    //     // DEBUG_PRINT("ADC disabled\n");
    //     DMA_ClearITPendingBit(DMA_Stream, DMA_IT_HTIF0);
    // }

    if (DMA_GetITStatus(DMA_Stream, DMA_IT_TCIF4)) //&& (xSemaphoreTake(semaphoreHalfBuffer, 0) == pdTRUE))
    {
        // ADC_Cmd(ADC_n, DISABLE);
        // DEBUG_PRINT("ADC disabled\n");
        DMA_ClearITPendingBit(DMA_Stream, DMA_IT_TCIF4);
        ADC_Done = 1;
    }
}

static void magneticInit()
{
    if (isInit)
    {
        return;
    }
    DEBUG_PRINT("MAGNETIC init started!\n");

    xTaskCreate(mytask, MAGNETIC_TASK_NAME,
                MAGNETIC_TASK_STACKSIZE, NULL, MAGNETIC_TASK_PRI, NULL);

    // DEBUG_PRINT("MAGNETIC init ended!\n");
    isInit = true;
}

static bool magneticTest()
{
    DEBUG_PRINT("MAGNETIC test passed!\n");
    return true;
}

static const DeckDriver magneticDriver = {
    .name = "Magnetic",
    .usedGpio = DECK_USING_PA3,
    .init = magneticInit,
    .test = magneticTest,
};

DECK_DRIVER(magneticDriver);
#define CONFIG_DEBUG_LOG_ENABLE = y
LOG_GROUP_START(example)
LOG_ADD_DEBUG(LOG_UINT16, firstValue, &firstValue)
LOG_ADD_DEBUG(LOG_UINT16, FirstVolt, &FirstVolt)
LOG_ADD_DEBUG(LOG_FLOAT, logNeroAmpl, &logNeroAmpl)
LOG_ADD_DEBUG(LOG_FLOAT, maxval, &maxval)
LOG_GROUP_STOP(example)

PARAM_GROUP_START(FFT)
PARAM_ADD(PARAM_UINT32, IDx, &IDx)
PARAM_ADD(PARAM_UINT32, valueToExclude, &valueToExclude)
PARAM_GROUP_STOP(FFT)