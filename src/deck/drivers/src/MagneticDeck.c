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
#include "system.h"
#include "param.h"
#include "i2cdev.h"

#define CONFIG_DEBUG = y

// System HZ
#define SYSTEM_HZ 10
#define SYSTEM_PERIOD_MS (1000 / SYSTEM_HZ)

// ADC DMA configuration
#define ARRAY_SIZE 2048
uint32_t DMA_Buffer[ARRAY_SIZE];
ADC_TypeDef *ADC_n = ADC1;
DMA_Stream_TypeDef *DMA_Stream = DMA2_Stream4;
uint32_t DMA_Channel = DMA_Channel_0;
IRQn_Type DMA_IRQ = DMA2_Stream4_IRQn;
// 2^12 WHERE 12 IS THE NUMBER OF BITS OF THE MCU ADC = 4096
#define ADC_LEVELS 4096
#define ADC_MAX_VOLTAGE 3.0f
#define PCLK2 84e6f
#define ADC_PRESCALER 6.0f
// 12 from bit and 15 from the register value 12+15 = 27
#define ADC_Full_Sampling_Time 27.0f
#define Fc_ADC (PCLK2 / ADC_PRESCALER / ADC_Full_Sampling_Time)

// #define DMA_IRQ DMA2_Stream0_IRQn
uint8_t ADC_Channel = ADC_Channel_3;
volatile uint8_t test = 0;
static bool isInit = false;

// ADC flag to check if the conversion is done
volatile uint8_t ADC_Done = 0;

// FFT parameters
#define FFT_SIZE ARRAY_SIZE
float32_t fft_input[FFT_SIZE];
float32_t fft_output[FFT_SIZE];
float32_t fft_magnitude[FFT_SIZE / 2];
arm_rfft_fast_instance_f32 fft_instance;
uint16_t fft_length = FFT_SIZE;
#define BIN_SIZE (int)(Fc_ADC / FFT_SIZE)

// DAC Reference Voltage Params
#define V_REF_CRAZYFLIE 3.0f
#define DAC_BIT 12
#define DAC_LEVELS (1 << DAC_BIT) // 4096
#define DAC_STEP (VREF / DAC_LEVELS)
#define DECK_DAC_I2C_ADDRESS 0x4C //
#define DAC_WRITE_LENGTH 2
#define DAC_ANALOG_HW_DELAY_AFTER_SET 10
#define V_DD 3.3f

// Potentiometer Params
#define R10 200.0f  // 200 Ohm
#define RW_2_7V 155 // ohm  155 OHM TYPICAL_DC_WIPER_RESISTANCE
#define RW_5_5V 100 // ohm

#define POTENTIOMETER_ADDR 0x2F // 0x94 WRITE 0x95 READ
#define POTENTIOMETER_BIT 7
// #define POTENTIOMETER_STEPS pow(2, POTENTIOMETER_BIT) // 7 bit --> 128
#define POTENTIOMETER_NUMBER_OF_STEPS (1 << POTENTIOMETER_BIT) - 1 // 127
#define POTENTIOMETER_FULL_SCALE_RAB 50E3                          // 50 kOhm
#define POTENTIOMETER_ANALOG_HW_DELAY_AFTER_SET 10
float32_t GainValue_Setted = 0;
float GainValue = 10;

// ------ Anchors Parameters -------
// Resonance Freqs Anchors in Hz
#define NeroResFreq 213e3
#define NeroIdx (int)(NeroResFreq / BIN_SIZE)
#define Nero_M -2.804
#define Nero_Q -2.635
float Nero_Position[] = {-0.50, -0.25, 0};
#define Nero_Id 0

#define GialloResFreq 203e3
#define GialloIdx (int)(GialloResFreq / BIN_SIZE)
#define Giallo_M -2.887
#define Giallo_Q -2.629
float Giallo_Position[] = {+0.50, -0.25, 0};
#define Giallo_Id 1

#define GrigioResFreq 193e3
#define GrigioIdx (int)(GrigioResFreq / BIN_SIZE)
#define Grigio_M -2.902
#define Grigio_Q -2.647
float Grigio_Position[] = {+0.50, +0.25, 0};
#define Grigio_Id 2

#define RossoResFreq 183e3
#define RossoIdx (int)(RossoResFreq / BIN_SIZE)
#define Rosso_M -2.950
#define Rosso_Q -2.640
float Rosso_Position[] = {-0.50, +0.25, 0};
#define Rosso_Id 3

#define Default_MagneticStandardDeviation 0.0001f
volatile float32_t MagneticStandardDeviation = Default_MagneticStandardDeviation;
// -------  Debug variables -------

float32_t NeroAmpl = 0;
float32_t GialloAmpl = 0;
float32_t GrigioAmpl = 0;
float32_t RossoAmpl = 0;

// Distances
volatile float32_t Nero_distance = 0;
volatile float32_t Giallo_distance = 0;
volatile float32_t Grigio_distance = 0;
volatile float32_t Rosso_distance = 0;

// FFT
volatile uint16_t bin_size = BIN_SIZE;
volatile float Fc = Fc_ADC;
volatile uint16_t fft_size = FFT_SIZE;
volatile uint16_t Nero_Idx = NeroIdx;
volatile uint16_t Giallo_Idx = GialloIdx;
volatile uint16_t Grigio_Idx = GrigioIdx;
volatile uint16_t Rosso_Idx = RossoIdx;

// -------------------- DAC Functions -------------------------------
uint16_t ValueforDAC_from_DesideredVol(float desidered_Voltage)
{
    // DAC7571
    // https://www.ti.com/lit/ds/symlink/dac7571.pdf?ts=1718008558964&ref_url=https%253A%252F%252Fwww.ti.com%252Fproduct%252FDAC7571
    float V_out = desidered_Voltage;
    float D_raw = (V_out * DAC_LEVELS) / V_DD;
    // Round the value to the nearest integer value.
    uint16_t D = (uint16_t)roundf(D_raw);
    // check if the number is in [0-4095] otherwise raise and error
    if (D < 0 || D > DAC_LEVELS - 1)
    {
        DEBUG_PRINT("Error: The value of the DAC is not in the range [0-4095]:%d\n", D);
        return 0;
    }

    return D;
}

// ----------------- Potentiometer function ------------------------------
uint16_t WiperResistanceValue_From_Interpolation(float Vdd)
{
    // https://ww1.microchip.com/downloads/aemDocuments/documents/OTH/ProductDocuments/DataSheets/22147a.pdf

    float R_WB = ((5.5f - Vdd) / (5.5f - 2.7f)) * RW_2_7V + ((Vdd - 2.7f) / (5.5f - 2.7f)) * RW_5_5V;

    uint16_t D = (uint16_t)roundf(R_WB);
    return R_WB;
}

float Potentiometer_Resistance_Value_from_Desidered_Gain(float desidered_Gain)
{
    // https://ww1.microchip.com/downloads/aemDocuments/documents/OTH/ProductDocuments/DataSheets/22147a.pdf
    // this resistance is R12 in the schematic and it correspond to R_wb in formulas 6-2
    return (desidered_Gain - 1) * R10;
}

uint8_t Potentiometer_Value_To_Set(float desidered_Gain)
{
    // equation  6-2 datasheet  DS22147A
    // https://ww1.microchip.com/downloads/aemDocuments/documents/OTH/ProductDocuments/DataSheets/22147a.pdf

    float R_WB_from_Desidered_Gain = Potentiometer_Resistance_Value_from_Desidered_Gain(desidered_Gain);
    // DEBUG_PRINT("R_WB_from_Desidered_Gain: %f\n", R_WB_from_Desidered_Gain);
    // Step resistance (RS) is the resistance from one tap setting to the next. Values in  [4000-6000] Ohm, typical 5000 Ohm
    float R_S = POTENTIOMETER_FULL_SCALE_RAB / POTENTIOMETER_NUMBER_OF_STEPS;
    // Compute the Raw value of the N starting from equation 6-2
    uint16_t TYPICAL_DC_WIPER_RESISTANCE = WiperResistanceValue_From_Interpolation(V_DD);
    float N_raw = (R_WB_from_Desidered_Gain - TYPICAL_DC_WIPER_RESISTANCE) / R_S;
    // Round the value to the nearest integer value.
    uint8_t N = (uint8_t)roundf(N_raw);
    // check if the number is in [0-127] otherwise raise and error
    if (N > POTENTIOMETER_NUMBER_OF_STEPS)
    {
        DEBUG_PRINT("Error: The value of the potentiometer is not in the range [0-127]:%d\n", N);
        return 0;
    }
    return N;
}

// ---------------- DMA Interrupt Handler ------------------------------
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
        // DMA_Cmd(DMA_Stream, DISABLE);
        ADC_Cmd(ADC_n, DISABLE);
        // DEBUG_PRINT("ADC disabled\n");
        ADC_ContinuousModeCmd(ADC_n, DISABLE);
        DMA_ClearITPendingBit(DMA_Stream, DMA_IT_TCIF4);

        ADC_Done = 1;
    }
    if (DMA_GetITStatus(DMA_Stream, DMA_IT_TEIF4))
    {
        DEBUG_PRINT("DMA_IT_TEIF4\n");
        DMA_ClearITPendingBit(DMA_Stream, DMA_IT_TEIF4);
    }
}

void ADC1_IRQHandler(void)
{
    DEBUG_PRINT("ADC1_IRQHandler\n");
    if (ADC_GetFlagStatus(ADC_n, ADC_FLAG_OVR))
    {
        ADC_ClearFlag(ADC_n, ADC_FLAG_OVR);
    }
}

// -------------------- Cycling Function ----------------------------
// fft and kalman update function
void performFFT(uint32_t *Input_buffer_pointer, float32_t *Output_buffer_pointer, float32_t flattopCorrectionFactor)
{
    arm_q31_to_float((q31_t *)Input_buffer_pointer, Output_buffer_pointer, FFT_SIZE);
    int p = 0;
    for (p = 0; p < FFT_SIZE; p++)
    {
        // Fattore di scala della conversione della funzione, guarda doc per capire
        Output_buffer_pointer[p] = Output_buffer_pointer[p] * 2147483648.0f;

        // ADCLevelsToVolt
        Output_buffer_pointer[p] = Output_buffer_pointer[p] * (ADC_MAX_VOLTAGE / ADC_LEVELS);
    }

    // TODO: 2.4 Voltages max voltages before saturation, so if you see 2.4V in the output, then change the gain. Also viceversa, to be tested the mimumum value

    // apply flattop window
    arm_mult_f32(Output_buffer_pointer, (float32_t *)flattop_2048_lut, Output_buffer_pointer, FFT_SIZE);

    // perform FFT
    arm_rfft_fast_f32(&fft_instance, Output_buffer_pointer, fft_output, 0);

    // normalizing the fft output
    int i = 0;
    for (i = 0; i < FFT_SIZE; i++)
    {
        fft_output[i] = fft_output[i] / FFT_SIZE;
    }

    // computing the magnitude of the fft
    arm_cmplx_mag_f32(fft_output, fft_magnitude, FFT_SIZE / 2);

    // Extract the max for each anchor
    // NOTE:*2 is because the fft of a sin is 2 impulsive delta of A/2 amplitude therefore to get the full amplitude of the signal i need to multiply by 2
    // NOTE: *flattopCorrectionFactor is the correction factor for the flattop window that i applied

    // Extract the maximum value and its index around NeroIdx
    uint32_t maxindex;
    arm_max_f32(&fft_magnitude[NeroIdx - 2], 4, &NeroAmpl, &maxindex);
    NeroAmpl = NeroAmpl * flattopCorrectionFactor * 2.0f;

    // Calculate the maximum value and its index around GialloIdx
    arm_max_f32(&fft_magnitude[GialloIdx - 1], 3, &GialloAmpl, &maxindex);
    GialloAmpl = GialloAmpl * flattopCorrectionFactor * 2.0f;

    // Calculate the maximum value and its index around GrigioIdx
    arm_max_f32(&fft_magnitude[GrigioIdx - 1], 3, &GrigioAmpl, &maxindex);
    GrigioAmpl = GrigioAmpl * flattopCorrectionFactor * 2.0f;

    // Calculate the maximum value and its index around RossoIdx
    arm_max_f32(&fft_magnitude[RossoIdx - 1], 3, &RossoAmpl, &maxindex);
    RossoAmpl = RossoAmpl * flattopCorrectionFactor * 2.0f;

    // TO IMPLEMENT!!!!
    // check if the ADC is saturating
    // if (NeroAmpl == 2.4 || GialloAmpl == 0 || GrigioAmpl == 0 || RossoAmpl == 0)
    // {
    //     DEBUG_PRINT("ADC is saturating\n");
    //     return;
    // }

    // --------------------------------- 2D MEASURMENT MODEL - DISTANCE COMPUTATION ---------------------------------
    // THIS MODEL IS NOT USED TO UPDATE THE KALMAN FILTER, JUST AS CONFIRMATION

    // compute the the distances from the amplitude of each anchor
    Nero_distance = powf(10, (log10(NeroAmpl) - Nero_Q) / Nero_M);
    Giallo_distance = powf(10, (log10(GialloAmpl) - Giallo_Q) / Giallo_M);
    Grigio_distance = powf(10, (log10(GrigioAmpl) - Grigio_Q) / Grigio_M);
    Rosso_distance = powf(10, (log10(RossoAmpl) - Rosso_Q) / Rosso_M);

    // Nero
    distanceMeasurement_t dist_Nero;
    dist_Nero.distance = Nero_distance;
    dist_Nero.x = Nero_Position[0];
    dist_Nero.y = Nero_Position[1];
    dist_Nero.z = Nero_Position[2];
    dist_Nero.anchorId = Nero_Id;
    dist_Nero.stdDev = MagneticStandardDeviation;
    // DEBUG_PRINT("Nero Distance: %f\n", Nero_distance);
    // estimatorEnqueueDistance(&dist_Nero);

    // // Giallo
    distanceMeasurement_t dist_Giallo;
    dist_Giallo.distance = Giallo_distance;
    dist_Giallo.x = Giallo_Position[0];
    dist_Giallo.y = Giallo_Position[1];
    dist_Giallo.z = Giallo_Position[2];
    dist_Giallo.anchorId = Giallo_Id;
    dist_Giallo.stdDev = MagneticStandardDeviation;
    // estimatorEnqueueDistance(&dist_Giallo);

    // // Grigio
    distanceMeasurement_t dist_Grigio;
    dist_Grigio.distance = Grigio_distance;
    dist_Grigio.x = Grigio_Position[0];
    dist_Grigio.y = Grigio_Position[1];
    dist_Grigio.z = Grigio_Position[2];
    dist_Grigio.anchorId = Grigio_Id;
    dist_Grigio.stdDev = MagneticStandardDeviation;
    // estimatorEnqueueDistance(&dist_Grigio);

    // // Rosso
    distanceMeasurement_t dist_Rosso;
    dist_Rosso.distance = Rosso_distance;
    dist_Rosso.x = Rosso_Position[0];
    dist_Rosso.y = Rosso_Position[1];
    dist_Rosso.z = Rosso_Position[2];
    dist_Rosso.anchorId = Rosso_Id;
    dist_Rosso.stdDev = MagneticStandardDeviation;
    // estimatorEnqueueDistance(&dist_Rosso);

    // ---------------3D MEASUREMENT MODEL - POSITION COMPUTATION------------------------------
    // THIS MODEL IS USED TO UPDATE THE KALMAN FILTER

    voltMeasurement_t volt;

    volt.x[0] = Nero_Position[0];
    volt.y[0] = Nero_Position[1];
    volt.z[0] = Nero_Position[2];
    volt.stdDev[0] = MagneticStandardDeviation;
    volt.measuredVolt[0] = NeroAmpl;
    volt.anchorId[0] = Nero_Id;
    volt.resonanceFrequency[0] = NeroResFreq;

    volt.x[1] = Giallo_Position[0];
    volt.y[1] = Giallo_Position[1];
    volt.z[1] = Giallo_Position[2];
    volt.stdDev[1] = MagneticStandardDeviation;
    volt.measuredVolt[1] = GialloAmpl;
    volt.anchorId[1] = Giallo_Id;
    volt.resonanceFrequency[1] = GialloResFreq;

    volt.x[2] = Grigio_Position[0];
    volt.y[2] = Grigio_Position[1];
    volt.z[2] = Grigio_Position[2];
    volt.stdDev[2] = MagneticStandardDeviation;
    volt.measuredVolt[2] = GrigioAmpl;
    volt.anchorId[2] = Grigio_Id;
    volt.resonanceFrequency[2] = GrigioResFreq;

    volt.x[3] = Rosso_Position[0];
    volt.y[3] = Rosso_Position[1];
    volt.z[3] = Rosso_Position[2];
    volt.stdDev[3] = MagneticStandardDeviation;
    volt.measuredVolt[3] = RossoAmpl;
    volt.anchorId[3] = Rosso_Id;
    volt.resonanceFrequency[3] = RossoResFreq;

    estimatorEnqueueVolt(&volt);
}

static void mytask(void *param)
{
    DEBUG_PRINT("Wait for system starting\n");
    systemWaitStart();
    DEBUG_PRINT("System Started\n");

    // DAC Setup
    // SETTING the reference voltage for the DAC using i2c
    uint16_t ValueforDAC = ValueforDAC_from_DesideredVol(V_REF_CRAZYFLIE / 2);
    uint8_t DAC_Write[2] = {0, 0};
    DAC_Write[1] = ValueforDAC & 0x00FF;
    DAC_Write[0] = (ValueforDAC >> 8) & 0x00FF;
    DEBUG_PRINT("Value for DAC: %d\n", ValueforDAC);
    uint8_t dacWriteResult = i2cdevWrite(I2C1_DEV, DECK_DAC_I2C_ADDRESS, DAC_WRITE_LENGTH, &DAC_Write[0]);
    // wait some time that the Hw is setted
    vTaskDelay(M2T(DAC_ANALOG_HW_DELAY_AFTER_SET));
    DEBUG_PRINT("DAC Write result: %d\n", dacWriteResult);

    // Potentiometer Setup
    uint8_t Potentiometer_Value = Potentiometer_Value_To_Set(GainValue);
    DEBUG_PRINT("Potentiometer Value to set: %d\n", Potentiometer_Value);
    uint8_t result_read_potentiometer = i2cdevWrite(I2C1_DEV, POTENTIOMETER_ADDR, 1, &Potentiometer_Value);
    vTaskDelay(M2T(POTENTIOMETER_ANALOG_HW_DELAY_AFTER_SET));
    DEBUG_PRINT("Potentiometer Write result: %d\n", result_read_potentiometer);
    GainValue_Setted = GainValue;
    GainValue = 0;

    // reading the value of the potentiometer
    uint8_t PotentiometerReadValue = 10;
    uint8_t result_read = i2cdevRead(I2C1_DEV, POTENTIOMETER_ADDR, 1, &PotentiometerReadValue);

    vTaskDelay(M2T(POTENTIOMETER_ANALOG_HW_DELAY_AFTER_SET));
    DEBUG_PRINT("Potentiometer Readed result: %d\n", result_read);
    DEBUG_PRINT("Potentiometer Read Value: %d\n", PotentiometerReadValue);

    // gpio init
    GPIO_init_analog(DECK_GPIO_RX2);
    // DMA init
    DMA_inititalization(RCC_AHB1Periph_DMA2, DMA_Stream, DMA_Buffer, ADC_n, DMA_Channel, DMA_IRQ, ARRAY_SIZE);
    // adc init
    ADC_init_DMA_mode(RCC_APB2Periph_ADC1, ADC_n);
    // Call the ADC_DMA_start function
    ADC_DMA_start(ADC_n, ADC_Channel, 1, ADC_SampleTime_15Cycles);
    // Initialize the FFT instance
    arm_rfft_fast_init_f32(&fft_instance, fft_length);

    // Flattop correction factor calculation
    float32_t sum = 0.0;
    for (int i = 0; i < ARRAY_SIZE; i++)
    {
        sum += flattop_2048_lut[i];
    }
    float32_t flattopCorrectionFactor = ARRAY_SIZE / sum;
    DEBUG_PRINT("Flattop Correction Factor: %f\n", flattopCorrectionFactor);

    while (1)
    {
        if (GainValue > 0)
        {
            DEBUG_PRINT("UPDATE FROM USER ON GAIN!!!!\n");
            // Potentiometer Setup
            uint8_t Potentiometer_Value = Potentiometer_Value_To_Set(GainValue);
            DEBUG_PRINT("Potentiometer Value to set: %d\n", Potentiometer_Value);
            uint8_t result_read_potentiometer = i2cdevWrite(I2C1_DEV, POTENTIOMETER_ADDR, 1, &Potentiometer_Value);
            vTaskDelay(M2T(POTENTIOMETER_ANALOG_HW_DELAY_AFTER_SET));
            DEBUG_PRINT("Potentiometer Write result: %d\n", result_read_potentiometer);
            GainValue_Setted = GainValue;
            GainValue = 0;
        }
        if (ADC_Done == 1 && GainValue == 0)
        {
            // The DMA buffer is full, perform the FFT

            performFFT(DMA_Buffer, fft_input, flattopCorrectionFactor);

            // The FFT is done, restart the ADC
            DMA_Cmd(DMA_Stream, ENABLE);
            ADC_DMA_start(ADC_n, ADC_Channel, 1, ADC_SampleTime_15Cycles);

            ADC_Done = 0;
        }
        else
        {
            // SOMETIME THE ADC IS BLOCKED, SO I NEED TO RESTART IT MANUALLY, SEEMS NOT APPENING ANYMORE
            // ------------DEBUGGING STUFF--------------

            // NON È NESSUNA DI QUESTE CONDIZIONI

            // if (ADC_GetFlagStatus(ADC1, ADC_FLAG_OVR) != RESET)
            // {
            //     // Overrun occurred
            //     // Take appropriate actions to handle the overrun
            //     DEBUG_PRINT("Overrun occurred\n");
            //     ADC_ClearFlag(ADC1, ADC_FLAG_OVR); // Clear the overrun flag
            //     ADC_Done = 1;
            // }
            //             ADC_ClearFlag(ADC_n, ADC_FLAG_AWD);
            // ADC_ClearFlag(ADC_n, ADC_FLAG_EOC);
            // ADC_ClearFlag(ADC_n, ADC_FLAG_JEOC);
            // ADC_ClearFlag(ADC_n, ADC_FLAG_JSTRT);
            // ADC_ClearFlag(ADC_n, ADC_FLAG_STRT);
            // ADC_ClearFlag(ADC_n, ADC_FLAG_OVR);
            // if (ADC_GetFlagStatus(ADC1, ADC_FLAG_STRT) != RESET)
            // {
            //     DEBUG_PRINT("ADC_FLAG_STRT\n");
            //     ADC_ClearFlag(ADC1, ADC_FLAG_STRT);
            //     ADC_Done = 1;
            // }
            // if (ADC_GetFlagStatus(ADC1, ADC_FLAG_EOC) != RESET)
            // {
            //     DEBUG_PRINT("ADC_FLAG_EOC\n");
            //     ADC_ClearFlag(ADC1, ADC_FLAG_EOC);
            //     ADC_Done = 1;
            // }
            // if (ADC_GetFlagStatus(ADC1, ADC_FLAG_JEOC) != RESET)
            // {
            //     DEBUG_PRINT("ADC_FLAG_JEOC\n");
            //     ADC_ClearFlag(ADC1, ADC_FLAG_JEOC);
            //     ADC_Done = 1;
            // }
            // if (ADC_GetFlagStatus(ADC1, ADC_FLAG_JSTRT) != RESET)
            // {
            //     DEBUG_PRINT("ADC_FLAG_JSTRT\n");
            //     ADC_ClearFlag(ADC1, ADC_FLAG_JSTRT);
            //     ADC_Done = 1;
            // }
            // if (ADC_GetFlagStatus(ADC1, ADC_FLAG_AWD) != RESET)
            // {
            //     DEBUG_PRINT("ADC_FLAG_AWD\n");
            //     ADC_ClearFlag(ADC1, ADC_FLAG_AWD);
            //     ADC_Done = 1;
            // }
            // if (ADC_GetFlagStatus(ADC1, ADC_FLAG_OVR) != RESET)
            // {
            //     DEBUG_PRINT("ADC_FLAG_OVR\n");
            //     ADC_ClearFlag(ADC1, ADC_FLAG_OVR);
            //     ADC_Done = 1;
            // }

            // Se si blocca l'adc si puo' usare questo
            // non ha senso!!!
            ADC_Done = 1;
            DEBUG_PRINT("ADC_Done: %d\n", ADC_Done);
            // ma funziona cosi, non so poi come siano i dati
        }

        // Defining the delay between the executions
        vTaskDelay(M2T(10));
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

    isInit = true;
}

static bool magneticTest()
{
    DEBUG_PRINT("MAGNETIC test passed!\n");
    return true;
}

static const DeckDriver magneticDriver = {
    .name = "MagneticDeck",
    .usedGpio = DECK_USING_PA3,
    .usedPeriph = DECK_USING_I2C,
    .init = magneticInit,
    .test = magneticTest,
};

DECK_DRIVER(magneticDriver);

#define CONFIG_DEBUG_LOG_ENABLE = y

LOG_GROUP_START(MAGNETIC_DISTANCES)
LOG_ADD(LOG_FLOAT, Nero, &Nero_distance)
LOG_ADD(LOG_FLOAT, Giallo, &Giallo_distance)
LOG_ADD(LOG_FLOAT, Grigio, &Grigio_distance)
LOG_ADD(LOG_FLOAT, Rosso, &Rosso_distance)
LOG_GROUP_STOP(MAGNETIC_DISTANCES)

LOG_GROUP_START(MAGNETIC_VOLTAGES)
LOG_ADD(LOG_FLOAT, Nero, &NeroAmpl)
LOG_ADD(LOG_FLOAT, Giallo, &GialloAmpl)
LOG_ADD(LOG_FLOAT, Grigio, &GrigioAmpl)
LOG_ADD(LOG_FLOAT, Rosso, &RossoAmpl)
LOG_GROUP_STOP(MAGNETIC_VOLTAGES)

PARAM_GROUP_START(MAGNETIC_Params)
// PARAM_ADD_CORE(PARAM_UINT16, NeroResFreq, &NeroResFreq)
// volatile int Nero_IDX = NeroIdx;
PARAM_ADD(PARAM_FLOAT, std_magn, &MagneticStandardDeviation)
// PARAM_ADD_CORE(PARAM_UINT16, Nero_Index, &Nero_IDX)
// PARAM_ADD_CORE(PARAM_UINT16, Giallo_Index, &Giallo_Idx)
// PARAM_ADD_CORE(PARAM_UINT16, Grigio_Index, &Grigio_Idx)
// PARAM_ADD_CORE(PARAM_UINT16, Rosso_Index, &Rosso_Idx)
PARAM_GROUP_STOP(MAGNETIC_Params)

//// POTENTIOMETER
LOG_GROUP_START(Potentiometer_g_L)
LOG_ADD(LOG_FLOAT, GpS, &GainValue_Setted)
LOG_GROUP_STOP(Potentiometer_G_L)

PARAM_GROUP_START(Potentiometer_G_P)
PARAM_ADD(PARAM_FLOAT, G_pot, &GainValue)
PARAM_GROUP_STOP(Potentiometer_G_P)