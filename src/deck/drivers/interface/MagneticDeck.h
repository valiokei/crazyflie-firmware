#ifndef __MAGNETIC_DECK_H__
#define __MAGNETIC_DECK_H__

#include "stm32f4xx_adc.h"
#include "stm32f4xx_dma.h"
#include "arm_math.h"
#include "deck_analog.h"
#include "stm32f4xx.h"
#include "stm32f4xx_gpio.h"

#define MODEL_TO_USE 0 // 0 for Nelder-Mead strictly, 1 for Nelder-Mead with Kalman Filter
#define Z_MEASUREMENTS_TO_ACCUMULATE 5

// -------------------------- ADC -------------------------------------------
// ADC DMA configuration
#define ARRAY_SIZE 2048

// 2^12 WHERE 12 IS THE NUMBER OF BITS OF THE MCU ADC = 4096
#define ADC_LEVELS 4096
#define ADC_MAX_VOLTAGE 3.0f
#define PCLK2 84e6f
#define ADC_PRESCALER 6.0f
// 12 from bit and 15 from the register value 12+15 = 27
#define ADC_Full_Sampling_Time 27.0f
#define Fc_ADC (PCLK2 / ADC_PRESCALER / ADC_Full_Sampling_Time)

#define ADC_Channel_Default ADC_Channel_3;
// -------------------------- DMA -------------------------------------------
#define DMA_IRQ DMA2_Stream4_IRQn
#define MY_DMA_Channel DMA_Channel_0
#define MY_DMA_Stream DMA2_Stream4

// IRQn_Type DMA_IRQ = DMA2_Stream4_IRQn;

// -------------------------- FFT -------------------------------------------
#define FFT_SIZE ARRAY_SIZE
#define BIN_SIZE (int)(Fc_ADC / FFT_SIZE)
#define SATURATION_TRESHOLD 1.2f

// ------------------------ Measurement Model Params -------------------------------------------
#define Default_MagneticStandardDeviation 0.0001f
#define G_INA 1000.0f
#define Optimization_Model_STD 0.08f

#define offsetCoil 0.025f // 3 cm

// ------------------------ Linear Kalman Filter Params -------------------------------------------
#define STATE_DIM 3   // Stato: [x, y, z]
#define MEASURE_DIM 3 // Misure: [x, y, z]
#define EPSILON 1e-6f // Termine di regolarizzazione

typedef struct
{
    float F[STATE_DIM][STATE_DIM];
    float H[MEASURE_DIM][STATE_DIM];
    float Q[STATE_DIM][STATE_DIM];
    float R[MEASURE_DIM][MEASURE_DIM];
    float P[STATE_DIM][STATE_DIM];
    float x_state[STATE_DIM];
    float z_meas[MEASURE_DIM];

} KalmanFilter;

#define q_kf_default 0.025f
#define sigma_x_kf_default 0.05f // will be squared in the R matrix

// const float q = 0.05f;
// const float sigma_x = 0.05f * 0.05f;

// -------------------------  adaptive std on Measured Voltage -------------------------------
#define UseAdaptiveSTD 0
#define window_size 25

// ------------------------ Calibration -------------------------------------------
#define CALIBRATION_TIC_VALUE 500.0f // Number of measurements to perform the calibration

// ------------------------System HZ-------------------------------------------
#define SYSTEM_HZ 30
#define SYSTEM_PERIOD_MS (1000 / SYSTEM_HZ)

// ------------------------ Anchors Parameters -------------------------------------------
#define NUM_ANCHORS 4

// Resonance Freqs Anchors in Hz
// #define NeroResFreq 213e3
#define NeroResFreq 210e3
#define NeroIdx (int)(NeroResFreq / BIN_SIZE)
#define Nero_M -2.804
#define Nero_Q -2.635
#define Nero_Position_x 0.48f
#define Nero_Position_y 0.48f
#define Nero_Position_z +0.30f - 0.11f
#define Nero_Id 0

// #define GialloResFreq 203e3
#define GialloResFreq 199e3
#define GialloIdx (int)(GialloResFreq / BIN_SIZE)
#define Giallo_M -2.887
#define Giallo_Q -2.629
#define Giallo_Position_x +0.48f
#define Giallo_Position_y -0.48f
#define Giallo_Position_z +0.30f - 0.11f
#define Giallo_Id 1

// #define GrigioResFreq 193e3
#define GrigioResFreq 189e3
#define GrigioIdx (int)(GrigioResFreq / BIN_SIZE)
#define Grigio_M -2.902
#define Grigio_Q -2.647
#define Grigio_Position_x -0.48f
#define Grigio_Position_y -0.48f
#define Grigio_Position_z +0.30f - 0.11f
#define Grigio_Id 2

// #define RossoResFreq 183e3
#define RossoResFreq 181e3
#define RossoIdx (int)(RossoResFreq / BIN_SIZE)
#define Rosso_M -2.950
#define Rosso_Q -2.640
#define Rosso_Position_x -0.48f
#define Rosso_Position_y +0.48f
#define Rosso_Position_z +0.30f - 0.11f
#define Rosso_Id 3

// ------------------------ PHYSICAL COIL -------------------------------------------
#define RAY 0.019f
#define N_WOUNDS 5.0f
#define COIL_SURFACE (RAY * RAY * PI)
#define CURRENT 0.5f // This maybe can be improved
#define MU_0 1.25663706212e-06f

// ------------------------ Gain -------------------------------------------
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
// in serie con l'INa c'e un operazionale che produce un guadagno di 10
#define DefaultPotentiometerValue G_INA / 10.0f
#define OpAmpGainValue 10.0f

// -------------------------- DAC -------------------------------------------
// DAC Reference Voltage Params
#define V_REF_CRAZYFLIE 3.0f
#define DAC_BIT 12
#define DAC_LEVELS (1 << DAC_BIT) // 4096
#define DAC_STEP (VREF / DAC_LEVELS)
#define DECK_DAC_I2C_ADDRESS 0x4C //
#define DAC_WRITE_LENGTH 2
#define DAC_ANALOG_HW_DELAY_AFTER_SET 10
#define V_DD 3.3f

// -------------------------- NELDER-MEAD structs -------------------------------------------
#define NM_OPTIMIZER_IMPLEMENTATION
#define NM_NO_DEBUG_LOG

// ========================== Function Definitions ==========================

// -------------------------- Magnetic Deck -------------------------------------------

uint16_t ValueforDAC_from_DesideredVolt(float desidered_Voltage);

uint16_t WiperResistanceValue_From_Interpolation(float Vdd);

float Potentiometer_Resistance_Value_from_Desidered_Gain(float desidered_Gain);

uint8_t Potentiometer_Value_To_Set(float desidered_Gain);

void DMA2_Stream4_IRQHandler(void);

void ADC1_IRQHandler(void);
typedef struct
{
    float NeroAmpl;
    float GialloAmpl;
    float GrigioAmpl;
    float RossoAmpl;
    float AllAmpl[NUM_ANCHORS];
} FFT_Amplitudes;

FFT_Amplitudes performFFT(uint32_t *Input_buffer_pointer, float32_t *Output_buffer_pointer, float32_t flattopCorrectionFactor);

void check_saturations(FFT_Amplitudes *amplitudes, int *Id_in_saturation, bool *there_is_saturation);

void finalizeCycle();

// -------------------------- Measurement Model  -------------------------------------------

void get_B_field_for_a_Anchor(float *anchor_pos,
                              const float *tag_pos,
                              float *tag_or_versor,
                              float *B_field);

float V_from_B(float *B_field, float *rx_versor, float resonanceFreq, float Gain);

void V_rx_derivate_x_function(float T_x, float T_y, float T_z, float W_x, float W_y, float W_z,
                              float N, float current, float S, float Gain,
                              float A0_x, float A0_y, float A0_z, float A1_x, float A1_y, float A1_z,
                              float A2_x, float A2_y, float A2_z, float A3_x, float A3_y, float A3_z,
                              float F_0, float F_1, float F_2, float F_3, float V_rx_derivative_x[4]);

void V_rx_derivate_y_function(float T_x, float T_y, float T_z, float W_x, float W_y, float W_z,
                              float N, float current, float S, float Gain,
                              float A0_x, float A0_y, float A0_z, float A1_x, float A1_y, float A1_z,
                              float A2_x, float A2_y, float A2_z, float A3_x, float A3_y, float A3_z,
                              float F_0, float F_1, float F_2, float F_3, float V_rx_derivate_y[4]);

void V_rx_derivate_z_function(float T_x, float T_y, float T_z, float W_x, float W_y, float W_z,
                              float N, float current, float S, float Gain,
                              float A0_x, float A0_y, float A0_z, float A1_x, float A1_y, float A1_z,
                              float A2_x, float A2_y, float A2_z, float A3_x, float A3_y, float A3_z,
                              float F_0, float F_1, float F_2, float F_3, float V_rx_derivate_z[4]);

// -------------------------- NELDER-MEAD Functions -------------------------------------------
const float myCostFunction_3A(int n, const float *x, void *arg);
const float myCostFunction_4A(int n, const float *x, void *arg);

// -------------------------- Math Utils -------------------------------------------

float dot_product(float *a, float *b, int length);

float euclidean_distance(float *a, float *b, int length);

void getversor(float *a, float *b, float *u, int length);

float computeSTD(float *data, int arrayDimension);
#define SQUARE(x) ((x) * (x))

#endif // __MAGNETIC_DECK_H__