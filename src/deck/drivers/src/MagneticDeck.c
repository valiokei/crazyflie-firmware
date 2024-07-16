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
#include "MagneticDeck.h"

#define CONFIG_DEBUG = y
static bool isInit = false;

// adc INITIALIZATION
// ADC flag to check if the conversion is done
uint8_t ADC_Done = 0;
static ADC_TypeDef *ADC_n = ADC1;
static uint8_t ADC_Channel = ADC_Channel_Default;

// DMA Initialization
static uint32_t DMA_Buffer[ARRAY_SIZE];
static DMA_Stream_TypeDef *DMA_Stream = MY_DMA_Stream;
static uint32_t DMA_Channel = MY_DMA_Channel;

// FFT parameters
static float32_t fft_input[FFT_SIZE];
static float32_t fft_output[FFT_SIZE];
static float32_t fft_magnitude[FFT_SIZE / 2];
static arm_rfft_fast_instance_f32 fft_instance;
static uint16_t fft_length = FFT_SIZE;

// Potentiometer Params
float32_t GainValue_Setted = 0;
float GainValue = DefaultPotentiometerValue;
// Total Gain
static float TotalGain = 0;

float MagneticStandardDeviation = Default_MagneticStandardDeviation;

// -------  Debug variables -------
static int counterSaturation = 0;

static float32_t NeroAmpl = 0;
static float32_t GialloAmpl = 0;
static float32_t GrigioAmpl = 0;
static float32_t RossoAmpl = 0;

// FFT
static uint16_t bin_size = BIN_SIZE;
static uint16_t fft_size = FFT_SIZE;
// static uint16_t Nero_Idx = NeroIdx;
// static uint16_t Giallo_Idx = GialloIdx;
// static uint16_t Grigio_Idx = GrigioIdx;
// static uint16_t Rosso_Idx = RossoIdx;

float skipThisMeasurement = 0;

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
    return D;
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

// --------------------------- Math Utils Functions ---------------------------
float dot_product(float *a, float *b, int length)
{
    float result = 0.0f;
    for (int i = 0; i < length; i++)
    {
        result += a[i] * b[i];
    }
    return result;
}

float euclidean_distance(float *a, float *b, int length)
{
    float sum = 0.0f;
    for (int i = 0; i < length; i++)
    {
        float diff = a[i] - b[i];
        sum += diff * diff;
    }
    return sqrtf(sum);
}

void getversor(float *a, float *b, float *u, int length)
{
    float d = euclidean_distance(a, b, length);
    for (int i = 0; i < length; i++)
    {
        u[i] = (a[i] - b[i]) / d;
    }
}

float computeSTD(float *data, int arrayDimension)
{
    float sum = 0.0;
    float mean = 0.0;
    float standardDeviation = 0.0;

    int i;

    for (i = 0; i < arrayDimension; i++)
    {
        sum += data[i];
    }

    mean = sum / arrayDimension;

    for (i = 0; i < arrayDimension; i++)
        standardDeviation += powf(data[i] - mean, 2);

    return sqrtf(standardDeviation / arrayDimension);
}

// ---------------- Measurement Model Functions ------------------------------
void get_B_field_for_a_Anchor(float *anchor_pos,
                              float *tag_pos,
                              float *tag_or_versor,
                              float *B_field)
{
    float tx_rx_versor[3];
    getversor(anchor_pos, tag_pos, tx_rx_versor, 3);
    // DEBUG_PRINT("tx_rx_versor = %f\n", tx_rx_versor);

    float magnetic_dipole_moment_tx[3] = {
        N_WOUNDS * COIL_SURFACE * CURRENT * tag_or_versor[0],
        N_WOUNDS * COIL_SURFACE * CURRENT * tag_or_versor[1],
        N_WOUNDS * COIL_SURFACE * CURRENT * tag_or_versor[2]};
    // DEBUG_PRINT("magnetic_dipole_moment_tx = %f\n", magnetic_dipole_moment_tx);

    float tx_rx_distance = euclidean_distance(anchor_pos, tag_pos, 3);

    // DEBUG_PRINT("tx_rx_distance = %f\n", tx_rx_distance);

    float dot_product_B_temp = dot_product(magnetic_dipole_moment_tx, tx_rx_versor, 3);

    float constant_Bfield_constant_1 = (MU_0 / (4.0f * PI)) / powf(tx_rx_distance, 3);
    float B_temp[3] = {
        constant_Bfield_constant_1 * (3.0f * dot_product_B_temp * tx_rx_versor[0] - magnetic_dipole_moment_tx[0]),
        constant_Bfield_constant_1 * (3.0f * dot_product_B_temp * tx_rx_versor[1] - magnetic_dipole_moment_tx[1]),
        constant_Bfield_constant_1 * (3.0f * dot_product_B_temp * tx_rx_versor[2] - magnetic_dipole_moment_tx[2])};

    // float magnetic_dipole_moment_tx_magnitude = euclidean_distance(magnetic_dipole_moment_tx, magnetic_dipole_moment_tx, 3);
    // compute the norm of the magnetic dipole moment
    // float magnetic_dipole_moment_tx_magnitude = sqrtf(magnetic_dipole_moment_tx[0] * magnetic_dipole_moment_tx[0] +
    //                                                   magnetic_dipole_moment_tx[1] * magnetic_dipole_moment_tx[1] +
    //                                                   magnetic_dipole_moment_tx[2] * magnetic_dipole_moment_tx[2]);

    // normalizing the B field (NOT NEEDED!!!!)
    B_field[0] = B_temp[0]; // * magnetic_dipole_moment_tx_magnitude;
    B_field[1] = B_temp[1]; // * magnetic_dipole_moment_tx_magnitude;
    B_field[2] = B_temp[2]; // * magnetic_dipole_moment_tx_magnitude;
}

float V_from_B(float *B_field, float *rx_versor, float resonanceFreq, float Gain)
{
    float dot_product_V = dot_product(B_field, rx_versor, 3);

    float V = Gain * fabsf(2.0f * PI * resonanceFreq * PI * RAY * RAY * N_WOUNDS * dot_product_V);
    return V;
}

void V_rx_derivate_x_function(float T_x, float T_y, float T_z, float W_x, float W_y, float W_z,
                              float N, float current, float S, float Gain,
                              float A0_x, float A0_y, float A0_z, float A1_x, float A1_y, float A1_z,
                              float A2_x, float A2_y, float A2_z, float A3_x, float A3_y, float A3_z,
                              float F_0, float F_1, float F_2, float F_3, float V_rx_derivative_x[4])
{
    float t2 = A0_x * 2.0f;
    float t3 = A1_x * 2.0f;
    float t4 = A2_x * 2.0f;
    float t5 = A3_x * 2.0f;
    float t6 = T_x * 2.0f;
    float t7 = 1.0f / PI;
    float t8 = -T_x;
    float t10 = -T_y;
    float t11 = -T_z;
    float t12 = N * S * W_x * current;
    float t13 = N * S * W_y * current;
    float t14 = N * S * W_z * current;
    float t9 = -t6;
    float t15 = A0_x + t8;
    float t16 = A1_x + t8;
    float t17 = A2_x + t8;
    float t18 = A3_x + t8;
    float t19 = A0_y + t10;
    float t20 = A1_y + t10;
    float t21 = A2_y + t10;
    float t22 = A3_y + t10;
    float t23 = A0_z + t11;
    float t24 = A1_z + t11;
    float t25 = A2_z + t11;
    float t26 = A3_z + t11;
    float t30 = t2 + t9;
    float t31 = t3 + t9;
    float t32 = t4 + t9;
    float t33 = t5 + t9;
    float t34 = t15 * t15;
    float t35 = t16 * t16;
    float t36 = t17 * t17;
    float t37 = t18 * t18;
    float t38 = t19 * t19;
    float t39 = t20 * t20;
    float t40 = t21 * t21;
    float t41 = t22 * t22;
    float t42 = t23 * t23;
    float t43 = t24 * t24;
    float t44 = t25 * t25;
    float t45 = t26 * t26;
    float t46 = t34 + t38 + t42;
    float t47 = t35 + t39 + t43;
    float t48 = t36 + t40 + t44;
    float t49 = t37 + t41 + t45;
    float t50 = 1.0f / sqrtf(t46);
    float t52 = 1.0f / sqrtf(t47);
    float t55 = 1.0f / sqrtf(t48);
    float t58 = 1.0f / sqrtf(t49);
    float t51 = t50 * t50 * t50;
    float t53 = t50 * t50 * t50 * t50 * t50;
    float t54 = t52 * t52 * t52;
    float t56 = t52 * t52 * t52 * t52 * t52;
    float t57 = t55 * t55 * t55;
    float t59 = t55 * t55 * t55 * t55 * t55;
    float t60 = t58 * t58 * t58;
    float t61 = t58 * t58 * t58 * t58 * t58;
    float t62 = t12 * t50 * 3.0f;
    float t63 = t12 * t52 * 3.0f;
    float t64 = t12 * t55 * 3.0f;
    float t65 = t12 * t58 * 3.0f;
    float t74 = t13 * t19 * t50 * 3.0f;
    float t75 = t13 * t20 * t52 * 3.0f;
    float t76 = t13 * t21 * t55 * 3.0f;
    float t77 = t13 * t22 * t58 * 3.0f;
    float t78 = t14 * t23 * t50 * 3.0f;
    float t79 = t14 * t24 * t52 * 3.0f;
    float t80 = t14 * t25 * t55 * 3.0f;
    float t81 = t14 * t26 * t58 * 3.0f;
    float t66 = -t62;
    float t67 = -t63;
    float t68 = -t64;
    float t69 = -t65;
    float t70 = t15 * t62;
    float t71 = t16 * t63;
    float t72 = t17 * t64;
    float t73 = t18 * t65;
    float t82 = t12 * t15 * t30 * t51 * (3.0f / 2.0f);
    float t83 = t12 * t16 * t31 * t54 * (3.0f / 2.0f);
    float t84 = t12 * t17 * t32 * t57 * (3.0f / 2.0f);
    float t85 = t12 * t18 * t33 * t60 * (3.0f / 2.0f);
    float t86 = t13 * t19 * t30 * t51 * (3.0f / 2.0f);
    float t87 = t13 * t20 * t31 * t54 * (3.0f / 2.0f);
    float t88 = t13 * t21 * t32 * t57 * (3.0f / 2.0f);
    float t89 = t13 * t22 * t33 * t60 * (3.0f / 2.0f);
    float t90 = t14 * t23 * t30 * t51 * (3.0f / 2.0f);
    float t91 = t14 * t24 * t31 * t54 * (3.0f / 2.0f);
    float t92 = t14 * t25 * t32 * t57 * (3.0f / 2.0f);
    float t93 = t14 * t26 * t33 * t60 * (3.0f / 2.0f);
    float t94 = t70 + t74 + t78;
    float t95 = t71 + t75 + t79;
    float t96 = t72 + t76 + t80;
    float t97 = t73 + t77 + t81;
    float t122 = t66 + t82 + t86 + t90;
    float t123 = t67 + t83 + t87 + t91;
    float t124 = t68 + t84 + t88 + t92;
    float t125 = t69 + t85 + t89 + t93;
    float t98 = t15 * t50 * t94;
    float t99 = t16 * t52 * t95;
    float t100 = t17 * t55 * t96;
    float t101 = t18 * t58 * t97;
    float t102 = t19 * t50 * t94;
    float t103 = t20 * t52 * t95;
    float t104 = t21 * t55 * t96;
    float t105 = t22 * t58 * t97;
    float t106 = t23 * t50 * t94;

    V_rx_derivative_x[0] = F_0 * Gain * N * S * PI *
                           copysignf(1.0f, MU_0 * W_x * t7 * t51 * (t12 - t98) * (-1.0f / 4.0f) -
                                               (MU_0 * W_y * t7 * t51 * (t13 - t102)) / 4.0f -
                                               (MU_0 * W_z * t7 * t51 * (t14 - t106)) / 4.0f) *
                           ((MU_0 * W_y * t7 * t51 * (t19 * t50 * t122 + (t19 * t30 * t51 * t94) / 2.0f)) / 4.0f +
                            (MU_0 * W_z * t7 * t51 * (t23 * t50 * t122 + (t23 * t30 * t51 * t94) / 2.0f)) / 4.0f +
                            (MU_0 * W_x * t7 * t51 * (-t50 * t94 + t15 * t50 * t122 + (t15 * t30 * t51 * t94) / 2.0f)) / 4.0f -
                            MU_0 * W_x * t7 * t30 * t53 * (t12 - t98) * (3.0f / 8.0f) -
                            MU_0 * W_y * t7 * t30 * t53 * (t13 - t102) * (3.0f / 8.0f) -
                            MU_0 * W_z * t7 * t30 * t53 * (t14 - t106) * (3.0f / 8.0f)) *
                           2.0f;

    float t107 = t24 * t52 * t95;
    V_rx_derivative_x[1] = F_1 * Gain * N * S * PI *
                           copysignf(1.0f, MU_0 * W_x * t7 * t54 * (t12 - t99) * (-1.0f / 4.0f) -
                                               (MU_0 * W_y * t7 * t54 * (t13 - t103)) / 4.0f -
                                               (MU_0 * W_z * t7 * t54 * (t14 - t107)) / 4.0f) *
                           ((MU_0 * W_y * t7 * t54 * (t20 * t52 * t123 + (t20 * t31 * t54 * t95) / 2.0f)) / 4.0f +
                            (MU_0 * W_z * t7 * t54 * (t24 * t52 * t123 + (t24 * t31 * t54 * t95) / 2.0f)) / 4.0f +
                            (MU_0 * W_x * t7 * t54 * (-t52 * t95 + t16 * t52 * t123 + (t16 * t31 * t54 * t95) / 2.0f)) / 4.0f -
                            MU_0 * W_x * t7 * t31 * t56 * (t12 - t99) * (3.0f / 8.0f) -
                            MU_0 * W_y * t7 * t31 * t56 * (t13 - t103) * (3.0f / 8.0f) -
                            MU_0 * W_z * t7 * t31 * t56 * (t14 - t107) * (3.0f / 8.0f)) *
                           2.0f;

    float t108 = t25 * t55 * t96;
    V_rx_derivative_x[2] = F_2 * Gain * N * S * PI *
                           copysignf(1.0f, MU_0 * W_x * t7 * t57 * (t12 - t100) * (-1.0f / 4.0f) -
                                               (MU_0 * W_y * t7 * t57 * (t13 - t104)) / 4.0f -
                                               (MU_0 * W_z * t7 * t57 * (t14 - t108)) / 4.0f) *
                           ((MU_0 * W_y * t7 * t57 * (t21 * t55 * t124 + (t21 * t32 * t57 * t96) / 2.0f)) / 4.0f +
                            (MU_0 * W_z * t7 * t57 * (t25 * t55 * t124 + (t25 * t32 * t57 * t96) / 2.0f)) / 4.0f +
                            (MU_0 * W_x * t7 * t57 * (-t55 * t96 + t17 * t55 * t124 + (t17 * t32 * t57 * t96) / 2.0f)) / 4.0f -
                            MU_0 * W_x * t7 * t32 * t59 * (t12 - t100) * (3.0f / 8.0f) -
                            MU_0 * W_y * t7 * t32 * t59 * (t13 - t104) * (3.0f / 8.0f) -
                            MU_0 * W_z * t7 * t32 * t59 * (t14 - t108) * (3.0f / 8.0f)) *
                           2.0f;

    float t109 = t26 * t58 * t97;
    V_rx_derivative_x[3] = F_3 * Gain * N * S * PI *
                           copysignf(1.0f, MU_0 * W_x * t7 * t60 * (t12 - t101) * (-1.0f / 4.0f) -
                                               (MU_0 * W_y * t7 * t60 * (t13 - t105)) / 4.0f -
                                               (MU_0 * W_z * t7 * t60 * (t14 - t109)) / 4.0f) *
                           ((MU_0 * W_y * t7 * t60 * (t22 * t58 * t125 + (t22 * t33 * t60 * t97) / 2.0f)) / 4.0f +
                            (MU_0 * W_z * t7 * t60 * (t26 * t58 * t125 + (t26 * t33 * t60 * t97) / 2.0f)) / 4.0f +
                            (MU_0 * W_x * t7 * t60 * (-t58 * t97 + t18 * t58 * t125 + (t18 * t33 * t60 * t97) / 2.0f)) / 4.0f -
                            MU_0 * W_x * t7 * t33 * t61 * (t12 - t101) * (3.0f / 8.0f) -
                            MU_0 * W_y * t7 * t33 * t61 * (t13 - t105) * (3.0f / 8.0f) -
                            MU_0 * W_z * t7 * t33 * t61 * (t14 - t109) * (3.0f / 8.0f)) *
                           2.0f;
}

void V_rx_derivate_y_function(float T_x, float T_y, float T_z, float W_x, float W_y, float W_z,
                              float N, float current, float S, float Gain,
                              float A0_x, float A0_y, float A0_z, float A1_x, float A1_y, float A1_z,
                              float A2_x, float A2_y, float A2_z, float A3_x, float A3_y, float A3_z,
                              float F_0, float F_1, float F_2, float F_3, float V_rx_derivate_y[4])
{
    float t2 = A0_y * 2.0f;
    float t3 = A1_y * 2.0f;
    float t4 = A2_y * 2.0f;
    float t5 = A3_y * 2.0f;
    float t6 = T_y * 2.0f;
    float t7 = 1.0f / PI;
    float t8 = -T_x;
    float t9 = -T_y;
    float t11 = -T_z;
    float t12 = N * S * W_x * current;
    float t13 = N * S * W_y * current;
    float t14 = N * S * W_z * current;
    float t10 = -t6;
    float t15 = A0_x + t8;
    float t16 = A1_x + t8;
    float t17 = A2_x + t8;
    float t18 = A3_x + t8;
    float t19 = A0_y + t9;
    float t20 = A1_y + t9;
    float t21 = A2_y + t9;
    float t22 = A3_y + t9;
    float t23 = A0_z + t11;
    float t24 = A1_z + t11;
    float t25 = A2_z + t11;
    float t26 = A3_z + t11;
    float t30 = t2 + t10;
    float t31 = t3 + t10;
    float t32 = t4 + t10;
    float t33 = t5 + t10;
    float t34 = t15 * t15;
    float t35 = t16 * t16;
    float t36 = t17 * t17;
    float t37 = t18 * t18;
    float t38 = t19 * t19;
    float t39 = t20 * t20;
    float t40 = t21 * t21;
    float t41 = t22 * t22;
    float t42 = t23 * t23;
    float t43 = t24 * t24;
    float t44 = t25 * t25;
    float t45 = t26 * t26;
    float t46 = t34 + t38 + t42;
    float t47 = t35 + t39 + t43;
    float t48 = t36 + t40 + t44;
    float t49 = t37 + t41 + t45;
    float t50 = 1.0f / sqrtf(t46);
    float t52 = 1.0f / sqrtf(t47);
    float t55 = 1.0f / sqrtf(t48);
    float t58 = 1.0f / sqrtf(t49);
    float t51 = t50 * t50 * t50;
    float t53 = t50 * t50 * t50 * t50 * t50;
    float t54 = t52 * t52 * t52;
    float t56 = t52 * t52 * t52 * t52 * t52;
    float t57 = t55 * t55 * t55;
    float t59 = t55 * t55 * t55 * t55 * t55;
    float t60 = t58 * t58 * t58;
    float t61 = t58 * t58 * t58 * t58 * t58;
    float t62 = t13 * t50 * 3.0f;
    float t63 = t13 * t52 * 3.0f;
    float t64 = t13 * t55 * 3.0f;
    float t65 = t13 * t58 * 3.0f;
    float t70 = t12 * t15 * t50 * 3.0f;
    float t71 = t12 * t16 * t52 * 3.0f;
    float t72 = t12 * t17 * t55 * 3.0f;
    float t73 = t12 * t18 * t58 * 3.0f;
    float t78 = t14 * t23 * t50 * 3.0f;
    float t79 = t14 * t24 * t52 * 3.0f;
    float t80 = t14 * t25 * t55 * 3.0f;
    float t81 = t14 * t26 * t58 * 3.0f;
    float t66 = -t62;
    float t67 = -t63;
    float t68 = -t64;
    float t69 = -t65;
    float t74 = t19 * t62;
    float t75 = t20 * t63;
    float t76 = t21 * t64;
    float t77 = t22 * t65;
    float t82 = t12 * t15 * t30 * t51 * (3.0f / 2.0f);
    float t83 = t12 * t16 * t31 * t54 * (3.0f / 2.0f);
    float t84 = t12 * t17 * t32 * t57 * (3.0f / 2.0f);
    float t85 = t12 * t18 * t33 * t60 * (3.0f / 2.0f);
    float t86 = t13 * t19 * t30 * t51 * (3.0f / 2.0f);
    float t87 = t13 * t20 * t31 * t54 * (3.0f / 2.0f);
    float t88 = t13 * t21 * t32 * t57 * (3.0f / 2.0f);
    float t89 = t13 * t22 * t33 * t60 * (3.0f / 2.0f);
    float t90 = t14 * t23 * t30 * t51 * (3.0f / 2.0f);
    float t91 = t14 * t24 * t31 * t54 * (3.0f / 2.0f);
    float t92 = t14 * t25 * t32 * t57 * (3.0f / 2.0f);
    float t93 = t14 * t26 * t33 * t60 * (3.0f / 2.0f);
    float t94 = t70 + t74 + t78;
    float t95 = t71 + t75 + t79;
    float t96 = t72 + t76 + t80;
    float t97 = t73 + t77 + t81;
    float t122 = t66 + t82 + t86 + t90;
    float t123 = t67 + t83 + t87 + t91;
    float t124 = t68 + t84 + t88 + t92;
    float t125 = t69 + t85 + t89 + t93;
    float t98 = t15 * t50 * t94;
    float t99 = t16 * t52 * t95;
    float t100 = t17 * t55 * t96;
    float t101 = t18 * t58 * t97;
    float t102 = t19 * t50 * t94;
    float t103 = t20 * t52 * t95;
    float t104 = t21 * t55 * t96;
    float t105 = t22 * t58 * t97;
    float t106 = t23 * t50 * t94;

    V_rx_derivate_y[0] = F_0 * Gain * N * S * PI *
                         copysignf(1.0f, MU_0 * W_x * t7 * t51 * (t12 - t98) * (-1.0f / 4.0f) -
                                             (MU_0 * W_y * t7 * t51 * (t13 - t102)) / 4.0f -
                                             (MU_0 * W_z * t7 * t51 * (t14 - t106)) / 4.0f) *
                         ((MU_0 * W_x * t7 * t51 * (t15 * t50 * t122 + (t15 * t30 * t51 * t94) / 2.0f)) / 4.0f +
                          (MU_0 * W_z * t7 * t51 * (t23 * t50 * t122 + (t23 * t30 * t51 * t94) / 2.0f)) / 4.0f +
                          (MU_0 * W_y * t7 * t51 * (-t50 * t94 + t19 * t50 * t122 + (t19 * t30 * t51 * t94) / 2.0f)) / 4.0f -
                          MU_0 * W_x * t7 * t30 * t53 * (t12 - t98) * (3.0f / 8.0f) -
                          MU_0 * W_y * t7 * t30 * t53 * (t13 - t102) * (3.0f / 8.0f) -
                          MU_0 * W_z * t7 * t30 * t53 * (t14 - t106) * (3.0f / 8.0f)) *
                         2.0f;

    float t107 = t24 * t52 * t95;
    V_rx_derivate_y[1] = F_1 * Gain * N * S * PI *
                         copysignf(1.0f, MU_0 * W_x * t7 * t54 * (t12 - t99) * (-1.0f / 4.0f) -
                                             (MU_0 * W_y * t7 * t54 * (t13 - t103)) / 4.0f -
                                             (MU_0 * W_z * t7 * t54 * (t14 - t107)) / 4.0f) *
                         ((MU_0 * W_x * t7 * t54 * (t16 * t52 * t123 + (t16 * t31 * t54 * t95) / 2.0f)) / 4.0f +
                          (MU_0 * W_z * t7 * t54 * (t24 * t52 * t123 + (t24 * t31 * t54 * t95) / 2.0f)) / 4.0f +
                          (MU_0 * W_y * t7 * t54 * (-t52 * t95 + t20 * t52 * t123 + (t20 * t31 * t54 * t95) / 2.0f)) / 4.0f -
                          MU_0 * W_x * t7 * t31 * t56 * (t12 - t99) * (3.0f / 8.0f) -
                          MU_0 * W_y * t7 * t31 * t56 * (t13 - t103) * (3.0f / 8.0f) -
                          MU_0 * W_z * t7 * t31 * t56 * (t14 - t107) * (3.0f / 8.0f)) *
                         2.0f;

    float t108 = t25 * t55 * t96;
    V_rx_derivate_y[2] = F_2 * Gain * N * S * PI *
                         copysignf(1.0f, MU_0 * W_x * t7 * t57 * (t12 - t100) * (-1.0f / 4.0f) -
                                             (MU_0 * W_y * t7 * t57 * (t13 - t104)) / 4.0f -
                                             (MU_0 * W_z * t7 * t57 * (t14 - t108)) / 4.0f) *
                         ((MU_0 * W_x * t7 * t57 * (t17 * t55 * t124 + (t17 * t32 * t57 * t96) / 2.0f)) / 4.0f +
                          (MU_0 * W_z * t7 * t57 * (t25 * t55 * t124 + (t25 * t32 * t57 * t96) / 2.0f)) / 4.0f +
                          (MU_0 * W_y * t7 * t57 * (-t55 * t96 + t21 * t55 * t124 + (t21 * t32 * t57 * t96) / 2.0f)) / 4.0f -
                          MU_0 * W_x * t7 * t32 * t59 * (t12 - t100) * (3.0f / 8.0f) -
                          MU_0 * W_y * t7 * t32 * t59 * (t13 - t104) * (3.0f / 8.0f) -
                          MU_0 * W_z * t7 * t32 * t59 * (t14 - t108) * (3.0f / 8.0f)) *
                         2.0f;

    float t109 = t26 * t58 * t97;
    V_rx_derivate_y[3] = F_3 * Gain * N * S * PI *
                         copysignf(1.0f, MU_0 * W_x * t7 * t60 * (t12 - t101) * (-1.0f / 4.0f) -
                                             (MU_0 * W_y * t7 * t60 * (t13 - t105)) / 4.0f -
                                             (MU_0 * W_z * t7 * t60 * (t14 - t109)) / 4.0f) *
                         ((MU_0 * W_x * t7 * t60 * (t18 * t58 * t125 + (t18 * t33 * t60 * t97) / 2.0f)) / 4.0f +
                          (MU_0 * W_z * t7 * t60 * (t26 * t58 * t125 + (t26 * t33 * t60 * t97) / 2.0f)) / 4.0f +
                          (MU_0 * W_y * t7 * t60 * (-t58 * t97 + t22 * t58 * t125 + (t22 * t33 * t60 * t97) / 2.0f)) / 4.0f -
                          MU_0 * W_x * t7 * t33 * t61 * (t12 - t101) * (3.0f / 8.0f) -
                          MU_0 * W_y * t7 * t33 * t61 * (t13 - t105) * (3.0f / 8.0f) -
                          MU_0 * W_z * t7 * t33 * t61 * (t14 - t109) * (3.0f / 8.0f)) *
                         2.0f;
}

void V_rx_derivate_z_function(float T_x, float T_y, float T_z, float W_x, float W_y, float W_z,
                              float N, float current, float S, float Gain,
                              float A0_x, float A0_y, float A0_z, float A1_x, float A1_y, float A1_z,
                              float A2_x, float A2_y, float A2_z, float A3_x, float A3_y, float A3_z,
                              float F_0, float F_1, float F_2, float F_3, float V_rx_derivate_z[4])
{
    float t2 = A0_z * 2.0f;
    float t3 = A1_z * 2.0f;
    float t4 = A2_z * 2.0f;
    float t5 = A3_z * 2.0f;
    float t6 = T_z * 2.0f;
    float t7 = 1.0f / PI;
    float t8 = -T_x;
    float t9 = -T_y;
    float t10 = -T_z;
    float t12 = N * S * W_x * current;
    float t13 = N * S * W_y * current;
    float t14 = N * S * W_z * current;
    float t11 = -t6;
    float t15 = A0_x + t8;
    float t16 = A1_x + t8;
    float t17 = A2_x + t8;
    float t18 = A3_x + t8;
    float t19 = A0_y + t9;
    float t20 = A1_y + t9;
    float t21 = A2_y + t9;
    float t22 = A3_y + t9;
    float t23 = A0_z + t10;
    float t24 = A1_z + t10;
    float t25 = A2_z + t10;
    float t26 = A3_z + t10;
    float t30 = t2 + t11;
    float t31 = t3 + t11;
    float t32 = t4 + t11;
    float t33 = t5 + t11;
    float t34 = t15 * t15;
    float t35 = t16 * t16;
    float t36 = t17 * t17;
    float t37 = t18 * t18;
    float t38 = t19 * t19;
    float t39 = t20 * t20;
    float t40 = t21 * t21;
    float t41 = t22 * t22;
    float t42 = t23 * t23;
    float t43 = t24 * t24;
    float t44 = t25 * t25;
    float t45 = t26 * t26;
    float t46 = t34 + t38 + t42;
    float t47 = t35 + t39 + t43;
    float t48 = t36 + t40 + t44;
    float t49 = t37 + t41 + t45;
    float t50 = 1.0f / sqrtf(t46);
    float t52 = 1.0f / sqrtf(t47);
    float t55 = 1.0f / sqrtf(t48);
    float t58 = 1.0f / sqrtf(t49);
    float t51 = t50 * t50 * t50;
    float t53 = t50 * t50 * t50 * t50 * t50;
    float t54 = t52 * t52 * t52;
    float t56 = t52 * t52 * t52 * t52 * t52;
    float t57 = t55 * t55 * t55;
    float t59 = t55 * t55 * t55 * t55 * t55;
    float t60 = t58 * t58 * t58;
    float t61 = t58 * t58 * t58 * t58 * t58;
    float t62 = t14 * t50 * 3.0f;
    float t63 = t14 * t52 * 3.0f;
    float t64 = t14 * t55 * 3.0f;
    float t65 = t14 * t58 * 3.0f;
    float t70 = t12 * t15 * t50 * 3.0f;
    float t71 = t12 * t16 * t52 * 3.0f;
    float t72 = t12 * t17 * t55 * 3.0f;
    float t73 = t12 * t18 * t58 * 3.0f;
    float t74 = t13 * t19 * t50 * 3.0f;
    float t75 = t13 * t20 * t52 * 3.0f;
    float t76 = t13 * t21 * t55 * 3.0f;
    float t77 = t13 * t22 * t58 * 3.0f;
    float t66 = -t62;
    float t67 = -t63;
    float t68 = -t64;
    float t69 = -t65;
    float t78 = t23 * t62;
    float t79 = t24 * t63;
    float t80 = t25 * t64;
    float t81 = t26 * t65;
    float t82 = t12 * t15 * t30 * t51 * (3.0f / 2.0f);
    float t83 = t12 * t16 * t31 * t54 * (3.0f / 2.0f);
    float t84 = t12 * t17 * t32 * t57 * (3.0f / 2.0f);
    float t85 = t12 * t18 * t33 * t60 * (3.0f / 2.0f);
    float t86 = t13 * t19 * t30 * t51 * (3.0f / 2.0f);
    float t87 = t13 * t20 * t31 * t54 * (3.0f / 2.0f);
    float t88 = t13 * t21 * t32 * t57 * (3.0f / 2.0f);
    float t89 = t13 * t22 * t33 * t60 * (3.0f / 2.0f);
    float t90 = t14 * t23 * t30 * t51 * (3.0f / 2.0f);
    float t91 = t14 * t24 * t31 * t54 * (3.0f / 2.0f);
    float t92 = t14 * t25 * t32 * t57 * (3.0f / 2.0f);
    float t93 = t14 * t26 * t33 * t60 * (3.0f / 2.0f);
    float t94 = t70 + t74 + t78;
    float t95 = t71 + t75 + t79;
    float t96 = t72 + t76 + t80;
    float t97 = t73 + t77 + t81;
    float t122 = t66 + t82 + t86 + t90;
    float t123 = t67 + t83 + t87 + t91;
    float t124 = t68 + t84 + t88 + t92;
    float t125 = t69 + t85 + t89 + t93;
    float t98 = t15 * t50 * t94;
    float t99 = t16 * t52 * t95;
    float t100 = t17 * t55 * t96;
    float t101 = t18 * t58 * t97;
    float t102 = t19 * t50 * t94;
    float t103 = t20 * t52 * t95;
    float t104 = t21 * t55 * t96;
    float t105 = t22 * t58 * t97;
    float t106 = t23 * t50 * t94;

    V_rx_derivate_z[0] = F_0 * Gain * N * S * PI *
                         copysignf(1.0f, MU_0 * W_x * t7 * t51 * (t12 - t98) * (-1.0f / 4.0f) -
                                             (MU_0 * W_y * t7 * t51 * (t13 - t102)) / 4.0f -
                                             (MU_0 * W_z * t7 * t51 * (t14 - t106)) / 4.0f) *
                         ((MU_0 * W_x * t7 * t51 * (t15 * t50 * t122 + (t15 * t30 * t51 * t94) / 2.0f)) / 4.0f +
                          (MU_0 * W_y * t7 * t51 * (t19 * t50 * t122 + (t19 * t30 * t51 * t94) / 2.0f)) / 4.0f +
                          (MU_0 * W_z * t7 * t51 * (-t50 * t94 + t23 * t50 * t122 + (t23 * t30 * t51 * t94) / 2.0f)) / 4.0f -
                          MU_0 * W_x * t7 * t30 * t53 * (t12 - t98) * (3.0f / 8.0f) -
                          MU_0 * W_y * t7 * t30 * t53 * (t13 - t102) * (3.0f / 8.0f) -
                          MU_0 * W_z * t7 * t30 * t53 * (t14 - t106) * (3.0f / 8.0f)) *
                         2.0f;

    float t107 = t24 * t52 * t95;
    V_rx_derivate_z[1] = F_1 * Gain * N * S * PI *
                         copysignf(1.0f, MU_0 * W_x * t7 * t54 * (t12 - t99) * (-1.0f / 4.0f) -
                                             (MU_0 * W_y * t7 * t54 * (t13 - t103)) / 4.0f -
                                             (MU_0 * W_z * t7 * t54 * (t14 - t107)) / 4.0f) *
                         ((MU_0 * W_x * t7 * t54 * (t16 * t52 * t123 + (t16 * t31 * t54 * t95) / 2.0f)) / 4.0f +
                          (MU_0 * W_y * t7 * t54 * (t20 * t52 * t123 + (t20 * t31 * t54 * t95) / 2.0f)) / 4.0f +
                          (MU_0 * W_z * t7 * t54 * (-t52 * t95 + t24 * t52 * t123 + (t24 * t31 * t54 * t95) / 2.0f)) / 4.0f -
                          MU_0 * W_x * t7 * t31 * t56 * (t12 - t99) * (3.0f / 8.0f) -
                          MU_0 * W_y * t7 * t31 * t56 * (t13 - t103) * (3.0f / 8.0f) -
                          MU_0 * W_z * t7 * t31 * t56 * (t14 - t107) * (3.0f / 8.0f)) *
                         2.0f;

    float t108 = t25 * t55 * t96;
    V_rx_derivate_z[2] = F_2 * Gain * N * S * PI *
                         copysignf(1.0f, MU_0 * W_x * t7 * t57 * (t12 - t100) * (-1.0f / 4.0f) -
                                             (MU_0 * W_y * t7 * t57 * (t13 - t104)) / 4.0f -
                                             (MU_0 * W_z * t7 * t57 * (t14 - t108)) / 4.0f) *
                         ((MU_0 * W_x * t7 * t57 * (t17 * t55 * t124 + (t17 * t32 * t57 * t96) / 2.0f)) / 4.0f +
                          (MU_0 * W_y * t7 * t57 * (t21 * t55 * t124 + (t21 * t32 * t57 * t96) / 2.0f)) / 4.0f +
                          (MU_0 * W_z * t7 * t57 * (-t55 * t96 + t25 * t55 * t124 + (t25 * t32 * t57 * t96) / 2.0f)) / 4.0f -
                          MU_0 * W_x * t7 * t32 * t59 * (t12 - t100) * (3.0f / 8.0f) -
                          MU_0 * W_y * t7 * t32 * t59 * (t13 - t104) * (3.0f / 8.0f) -
                          MU_0 * W_z * t7 * t32 * t59 * (t14 - t108) * (3.0f / 8.0f)) *
                         2.0f;

    float t109 = t26 * t58 * t97;
    V_rx_derivate_z[3] = F_3 * Gain * N * S * PI *
                         copysignf(1.0f, MU_0 * W_x * t7 * t60 * (t12 - t101) * (-1.0f / 4.0f) -
                                             (MU_0 * W_y * t7 * t60 * (t13 - t105)) / 4.0f -
                                             (MU_0 * W_z * t7 * t60 * (t14 - t109)) / 4.0f) *
                         ((MU_0 * W_x * t7 * t60 * (t18 * t58 * t125 + (t18 * t33 * t60 * t97) / 2.0f)) / 4.0f +
                          (MU_0 * W_y * t7 * t60 * (t22 * t58 * t125 + (t22 * t33 * t60 * t97) / 2.0f)) / 4.0f +
                          (MU_0 * W_z * t7 * t60 * (-t58 * t97 + t26 * t58 * t125 + (t26 * t33 * t60 * t97) / 2.0f)) / 4.0f -
                          MU_0 * W_x * t7 * t33 * t61 * (t12 - t101) * (3.0f / 8.0f) -
                          MU_0 * W_y * t7 * t33 * t61 * (t13 - t105) * (3.0f / 8.0f) -
                          MU_0 * W_z * t7 * t33 * t61 * (t14 - t109) * (3.0f / 8.0f)) *
                         2.0f;
}

// -------------------- Cycling Function ----------------------------

// fft and kalman update function
void performFFT(uint32_t *Input_buffer_pointer, float32_t *Output_buffer_pointer, float32_t flattopCorrectionFactor)
{
    if (skipThisMeasurement == 1)
    {
        return;
    }
    arm_q31_to_float((q31_t *)Input_buffer_pointer, Output_buffer_pointer, FFT_SIZE);
    int p = 0;
    for (p = 0; p < FFT_SIZE; p++)
    {
        // Fattore di scala della conversione della funzione, guarda doc per capire
        Output_buffer_pointer[p] = Output_buffer_pointer[p] * 2147483648.0f;

        // ADCLevelsToVolt
        Output_buffer_pointer[p] = Output_buffer_pointer[p] * (ADC_MAX_VOLTAGE / ADC_LEVELS);
    }

    // fin the maximum of the Output_buffer_pointer
    float maxSat;
    float maxIndexSat;
    arm_max_f32(Output_buffer_pointer, FFT_SIZE, &maxSat, &maxIndexSat);

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

    // --------------------------------- 2D MEASUREMENT MODEL - DISTANCE COMPUTATION ---------------------------------
    /*

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
    estimatorEnqueueDistance(&dist_Rosso);
    */

    // ---------------3D MEASUREMENT MODEL - POSITION COMPUTATION------------------------------
    // THIS MODEL IS USED TO UPDATE THE KALMAN FILTER

    voltMeasurement_t volt;

    volt.Id_in_saturation = 10;

    volt.x[0] = Nero_Position_x;
    volt.y[0] = Nero_Position_y;
    volt.z[0] = Nero_Position_z;
    volt.stdDev[0] = MagneticStandardDeviation;
    volt.measuredVolt[0] = NeroAmpl;
    volt.anchorId[0] = Nero_Id;
    volt.resonanceFrequency[0] = NeroResFreq;
    volt.GainValue = TotalGain;

    volt.x[1] = Giallo_Position_x;
    volt.y[1] = Giallo_Position_y;
    volt.z[1] = Giallo_Position_z;
    volt.stdDev[1] = MagneticStandardDeviation;
    volt.measuredVolt[1] = GialloAmpl;
    volt.anchorId[1] = Giallo_Id;
    volt.resonanceFrequency[1] = GialloResFreq;
    volt.GainValue = TotalGain;

    volt.x[2] = Grigio_Position_x;
    volt.y[2] = Grigio_Position_y;
    volt.z[2] = Grigio_Position_z;
    volt.stdDev[2] = MagneticStandardDeviation;
    volt.measuredVolt[2] = GrigioAmpl;
    volt.anchorId[2] = Grigio_Id;
    volt.resonanceFrequency[2] = GrigioResFreq;
    volt.GainValue = TotalGain;

    volt.x[3] = Rosso_Position_x;
    volt.y[3] = Rosso_Position_y;
    volt.z[3] = Rosso_Position_z;
    volt.stdDev[3] = MagneticStandardDeviation;
    volt.measuredVolt[3] = RossoAmpl;
    volt.anchorId[3] = Rosso_Id;
    volt.resonanceFrequency[3] = RossoResFreq;
    volt.GainValue = TotalGain;

    // check if the ADC is saturating
    // TODO: 2.4 Voltages max voltages before saturation, so if you see 2.4V in the output, then change the gain. Also viceversa, to be tested the mimumum value

    // ---------------------- SATURATION CASE -------------------------------------

    // check if the max is saturating
    if (maxSat >= 2.5f)
    {
        // if the max is saturating, then skip this measurement

        // if the max is saturating, then skip this measurement
        if (counterSaturation % 10 == 0)
        {
            DEBUG_PRINT("Saturating\n");
        }
        counterSaturation++;

        // finde the ancor more close to the max by looking at the max fft value
        // Find the maximum value among NeroAmpl, GialloAmpl, GrigioAmpl, RossoAmpl

        float maxAmpl = NeroAmpl;
        int maxAnchorId = Nero_Id;
        if (GialloAmpl > maxAmpl)
        {
            maxAmpl = GialloAmpl;
            maxAnchorId = Giallo_Id;
        }
        if (GrigioAmpl > maxAmpl)
        {
            maxAmpl = GrigioAmpl;
            maxAnchorId = Grigio_Id;
        }
        if (RossoAmpl > maxAmpl)
        {
            maxAmpl = RossoAmpl;
            maxAnchorId = Rosso_Id;
        }

        char *anchorName;
        switch (maxAnchorId)
        {
        case Nero_Id:
            anchorName = "Nero";
            break;
        case Giallo_Id:
            anchorName = "Giallo";
            break;
        case Grigio_Id:
            anchorName = "Grigio";
            break;
        case Rosso_Id:
            anchorName = "Rosso";
            break;
        default:
            anchorName = "Unknown";
            break;
        }
        // DEBUG_PRINT("Max Anchor Name: %s\n", anchorName);
        // Incresing the std for that measurement

        // volt.stdDev[maxAnchorId] = MagneticStandardDeviation * 2;

        volt.Id_in_saturation = maxAnchorId;
        // return;

        // DEBUG_PRINT("STD increased for %s *2\n", anchorName);
    }

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
    TotalGain = GainValue_Setted * OpAmpGainValue;

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
    DEBUG_PRINT("Flattop Correction Factor: %f\n", (double)flattopCorrectionFactor);

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
            TotalGain = GainValue_Setted * OpAmpGainValue;
        }
        if (ADC_Done == 1 && GainValue == 0)
        {
            // The DMA buffer is full, perform the FFT

            performFFT(DMA_Buffer, fft_input, flattopCorrectionFactor);

            // The FFT is done, restart the ADC
            DMA_Cmd(DMA_Stream, ENABLE);
            ADC_DMA_start(ADC_n, ADC_Channel, 1, ADC_SampleTime_15Cycles);
            skipThisMeasurement = 0;
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
            skipThisMeasurement = 1;
            // ma funziona cosi, non so poi come siano i dati
        }

        // Defining the delay between the executions
        vTaskDelay(M2T(SYSTEM_PERIOD_MS));
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

// LOG_GROUP_START(MAGNETIC_DISTANCES)
// LOG_ADD(LOG_FLOAT, Nero, &Nero_distance)
// LOG_ADD(LOG_FLOAT, Giallo, &Giallo_distance)
// LOG_ADD(LOG_FLOAT, Grigio, &Grigio_distance)
// LOG_ADD(LOG_FLOAT, Rosso, &Rosso_distance)
// LOG_GROUP_STOP(MAGNETIC_DISTANCES)

// LOG_GROUP_START(MAGNETIC_VOLTAGES)
// LOG_ADD(LOG_FLOAT, Nero, &NeroAmpl)
// LOG_ADD(LOG_FLOAT, Giallo, &GialloAmpl)
// LOG_ADD(LOG_FLOAT, Grigio, &GrigioAmpl)
// LOG_ADD(LOG_FLOAT, Rosso, &RossoAmpl)
// LOG_GROUP_STOP(MAGNETIC_VOLTAGES)

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
// LOG_GROUP_START(Potentiometer_g_L)
// LOG_ADD(LOG_FLOAT, GpS, &GainValue_Setted)
// LOG_GROUP_STOP(Potentiometer_G_L)

PARAM_GROUP_START(Potentiometer_G_P)
PARAM_ADD(PARAM_FLOAT, G_pot, &GainValue)
PARAM_GROUP_STOP(Potentiometer_G_P)