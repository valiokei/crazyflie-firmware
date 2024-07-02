#include "math.h"
#include "mm_magnetic_distance.h"
#include "log.h"
#include "debug.h"
#include "param.h"
#include "stabilizer_types.h"
#include "estimator.h"
#include "estimator_kalman.h"
#include <stdlib.h>

#define G_INA 100.0f
#define RAY 0.019f
#define N_WOUNDS 5.0f
#define COIL_SURFACE (RAY * RAY * PI)
#define CURRENT 0.5f
#define MU_0 1.25663706212e-06f
#define PI 3.14159265359f

// #define start_timer() *((volatile uint32_t *)0xE0001000) = 0x40000001 // Enable CYCCNT register
// #define stop_timer() *((volatile uint32_t *)0xE0001000) = 0x40000000  // Disable CYCCNT register
// #define get_timer() *((volatile uint32_t *)0xE0001004)

// logVarId_t gainID =  logGetVarId("Potentiometer_G_P", "GainValue");

// #define G_INA logGetFloat(GainID)

// uint32_t it1, it2; // start and stop flag
float MeasuredVoltages_calibrated[4] = {0.0f, 0.0f, 0.0f, 0.0f};
float MeasuredVoltages_raw[4] = {0.0f, 0.0f, 0.0f, 0.0f};
float PredictedVoltages[4] = {0.0f, 0.0f, 0.0f, 0.0f};
float Derivative_x[4] = {0.0f, 0.0f, 0.0f, 0.0f};
float Derivative_y[4] = {0.0f, 0.0f, 0.0f, 0.0f};
float Derivative_z[4] = {0.0f, 0.0f, 0.0f, 0.0f};
float B[4][3] = {{0.0f, 0.0f, 0.0f},
                 {0.0f, 0.0f, 0.0f},
                 {0.0f, 0.0f, 0.0f},
                 {0.0f, 0.0f, 0.0f}};
float E_1 = 0.0f;
float E_2 = 0.0f;
float E_3 = 0.0f;
float E_4 = 0.0f;

float T_pre_x = 0.0f;
float T_pre_y = 0.0f;
float T_pre_z = 0.0f;

float T_orientation[3] = {0.0f, 0.0f, 0.0f};

volatile float default_CG_a1 = 1.3f; //  NERO
volatile float CG_a1 = 1.0f;         //  NERO

volatile float default_CG_a2 = 0.635f; // GIALLO
volatile float CG_a2 = 1.0f;           // GIALLO

volatile float default_CG_a3 = 3.5f; // GRIGIO
volatile float CG_a3 = 1.0f;         // GRIGIO

volatile float default_CG_a4 = 7.5f; // ROSSO
volatile float CG_a4 = 1.0f;         // ROSSO

float CalibrationRapport_anchor1 = 1.0f;
float CalibrationRapport_anchor2 = 1.0f;
float CalibrationRapport_anchor3 = 1.0f;
float CalibrationRapport_anchor4 = 1.0f;

volatile uint32_t pos_idx = 0;

#define CALIBRATION_TIC_VALUE 100.0f // Replace with the desired constant value
volatile uint32_t currentCalibrationTick = 0;
static float calibrationMean[4][4] = {};

// DEBUG FUNCTION
double generate_gaussian_noise(double mu, double sigma)
{
    static const double epsilon = 1e-10;
    static const double two_pi = 2.0 * 3.14159265358979323846;

    static double z1;
    static int generate;
    generate = !generate;

    if (!generate)
        return z1 * sigma + mu;

    double u1, u2;
    do
    {
        u1 = rand() / (RAND_MAX + 1.0);
        u2 = rand() / (RAND_MAX + 1.0);
    } while (u1 <= epsilon);

    double z0;
    z0 = sqrt(-2.0 * log(u1)) * cos(two_pi * u2);
    z1 = sqrt(-2.0 * log(u1)) * sin(two_pi * u2);
    return z0 * sigma + mu;
}

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

void V_rx_derivate_x(float T_x, float T_y, float T_z, float W_x, float W_y, float W_z,
                     float N, float current, float S, float G_ina,
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

    V_rx_derivative_x[0] = F_0 * G_ina * N * S * PI *
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
    V_rx_derivative_x[1] = F_1 * G_ina * N * S * PI *
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
    V_rx_derivative_x[2] = F_2 * G_ina * N * S * PI *
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
    V_rx_derivative_x[3] = F_3 * G_ina * N * S * PI *
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

void V_rx_derivate_y(float T_x, float T_y, float T_z, float W_x, float W_y, float W_z,
                     float N, float current, float S, float G_ina,
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

    V_rx_derivate_y[0] = F_0 * G_ina * N * S * PI *
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
    V_rx_derivate_y[1] = F_1 * G_ina * N * S * PI *
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
    V_rx_derivate_y[2] = F_2 * G_ina * N * S * PI *
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
    V_rx_derivate_y[3] = F_3 * G_ina * N * S * PI *
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

void V_rx_derivate_z(float T_x, float T_y, float T_z, float W_x, float W_y, float W_z,
                     float N, float current, float S, float G_ina,
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

    V_rx_derivate_z[0] = F_0 * G_ina * N * S * PI *
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
    V_rx_derivate_z[1] = F_1 * G_ina * N * S * PI *
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
    V_rx_derivate_z[2] = F_2 * G_ina * N * S * PI *
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
    V_rx_derivate_z[3] = F_3 * G_ina * N * S * PI *
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

float V_from_B(float *B_field, float *rx_versor, float resonanceFreq)
{
    float dot_product_V = dot_product(B_field, rx_versor, 3);

    float V = G_INA * fabsf(2.0f * PI * resonanceFreq * PI * RAY * RAY * N_WOUNDS * dot_product_V);
    return V;
}

void kalmanCoreUpdateWithVolt(kalmanCoreData_t *this, voltMeasurement_t *voltAnchor)
{
    // DEBUG_PRINT("kalmanCoreUpdateWithVolt\n");

    float B_field_vector_1[3];
    float B_field_vector_2[3];
    float B_field_vector_3[3];
    float B_field_vector_4[3];

    if (currentCalibrationTick < CALIBRATION_TIC_VALUE)
    {
        calibrationMean[0][0] += voltAnchor->measuredVolt[0];
        // DEBUG_PRINT("measuredVolt[0] = %f\n", voltAnchor->measuredVolt[0]);

        calibrationMean[1][1] += voltAnchor->measuredVolt[1];
        // DEBUG_PRINT("measuredVolt[1] = %f\n", voltAnchor->measuredVolt[1]);

        calibrationMean[2][2] += voltAnchor->measuredVolt[2];
        // DEBUG_PRINT("measuredVolt[2] = %f\n", voltAnchor->measuredVolt[2]);

        calibrationMean[3][3] += voltAnchor->measuredVolt[3];
        // DEBUG_PRINT("measuredVolt[3] = %f\n", voltAnchor->measuredVolt[3]);

        currentCalibrationTick = currentCalibrationTick + 1;
        // DEBUG_PRINT("currentCalibrationTick = %d\n", currentCalibrationTick);
    }
    else
    {
        if (currentCalibrationTick == CALIBRATION_TIC_VALUE)
        {
            currentCalibrationTick = currentCalibrationTick + 1;

            float meanData_a1 = calibrationMean[0][0] / CALIBRATION_TIC_VALUE;
            float meanData_a2 = calibrationMean[1][1] / CALIBRATION_TIC_VALUE;
            float meanData_a3 = calibrationMean[2][2] / CALIBRATION_TIC_VALUE;
            float meanData_a4 = calibrationMean[3][3] / CALIBRATION_TIC_VALUE;

            DEBUG_PRINT("calibration data acquired!!\n");

            DEBUG_PRINT("calibration_Mean[0][0] = %f\n", (double)meanData_a1);
            DEBUG_PRINT("calibration_Mean[1][1] = %f\n", (double)meanData_a2);
            DEBUG_PRINT("calibration_Mean[2][2] = %f\n", (double)meanData_a3);
            DEBUG_PRINT("calibration_Mean[3][3] = %f\n", (double)meanData_a4);

            // computing the predicted voltages to calculate the rapport between the measured and the predicted
            float tag_pos_predicted_calibrated[3] = {this->S[KC_STATE_X], this->S[KC_STATE_Y], this->S[KC_STATE_Z]};
            T_pre_x = tag_pos_predicted_calibrated[0];
            T_pre_y = tag_pos_predicted_calibrated[1];
            T_pre_z = tag_pos_predicted_calibrated[2];

            // sarebbe il prodotto tra la matrice di rotazione e il versore  [0,0,1] iniziale
            // NON CAMBIA MAI!@!!!!!!!
            // float tag_or_versor_calibrated[3] = {this->R[0][2], this->R[1][2], this->R[2][2]};

            float RotationMatrix[3][3];
            // DEBUG_PRINT("RotationMatrix[0][0] = %f\n", (double)this->R[0][0]);
            // DEBUG_PRINT("RotationMatrix[0][1] = %f\n", (double)this->R[0][1]);
            // DEBUG_PRINT("RotationMatrix[0][2] = %f\n", (double)this->R[0][2]);
            // DEBUG_PRINT("RotationMatrix[1][0] = %f\n", (double)this->R[1][0]);
            // DEBUG_PRINT("RotationMatrix[1][1] = %f\n", (double)this->R[1][1]);
            // DEBUG_PRINT("RotationMatrix[1][2] = %f\n", (double)this->R[1][2]);
            // DEBUG_PRINT("RotationMatrix[2][0] = %f\n", RotationMatrix[2][0]);
            // DEBUG_PRINT("RotationMatrix[2][1] = %f\n", RotationMatrix[2][1]);
            // DEBUG_PRINT("RotationMatrix[2][2] = %f\n", RotationMatrix[2][2]);

            estimatorKalmanGetEstimatedRot((float *)RotationMatrix);
            float tag_or_versor_calibrated[3] = {RotationMatrix[0][2], RotationMatrix[1][2], RotationMatrix[2][2]};
            // logVarId_t idYaw = logGetVarId("stateEstimate", "yaw");
            // float yaw = 0.0f;
            // yaw = logGetFloat(idYaw);

            T_orientation[0] = tag_or_versor_calibrated[0];
            T_orientation[1] = tag_or_versor_calibrated[1];
            T_orientation[2] = tag_or_versor_calibrated[2];

            float anchor_1_pose[3] = {voltAnchor->x[0], voltAnchor->y[0], voltAnchor->z[0]};
            float anchor_2_pose[3] = {voltAnchor->x[1], voltAnchor->y[1], voltAnchor->z[1]};
            float anchor_3_pose[3] = {voltAnchor->x[2], voltAnchor->y[2], voltAnchor->z[2]};
            float anchor_4_pose[3] = {voltAnchor->x[3], voltAnchor->y[3], voltAnchor->z[3]};
            get_B_field_for_a_Anchor(anchor_1_pose, tag_pos_predicted_calibrated, tag_or_versor_calibrated, B_field_vector_1);
            get_B_field_for_a_Anchor(anchor_2_pose, tag_pos_predicted_calibrated, tag_or_versor_calibrated, B_field_vector_2);
            get_B_field_for_a_Anchor(anchor_3_pose, tag_pos_predicted_calibrated, tag_or_versor_calibrated, B_field_vector_3);
            get_B_field_for_a_Anchor(anchor_4_pose, tag_pos_predicted_calibrated, tag_or_versor_calibrated, B_field_vector_4);
            B[0][0] = B_field_vector_1[0];
            B[0][1] = B_field_vector_1[1];
            B[0][2] = B_field_vector_1[2];
            B[1][0] = B_field_vector_2[0];
            B[1][1] = B_field_vector_2[1];
            B[1][2] = B_field_vector_2[2];
            B[2][0] = B_field_vector_3[0];
            B[2][1] = B_field_vector_3[1];
            B[2][2] = B_field_vector_3[2];
            B[3][0] = B_field_vector_4[0];
            B[3][1] = B_field_vector_4[1];
            B[3][2] = B_field_vector_4[2];
            // computing the V_rx for each of the 4 anchors
            float V_rx_1 = V_from_B(B_field_vector_1, tag_or_versor_calibrated, voltAnchor->resonanceFrequency[0]);
            float V_rx_2 = V_from_B(B_field_vector_2, tag_or_versor_calibrated, voltAnchor->resonanceFrequency[1]);
            float V_rx_3 = V_from_B(B_field_vector_3, tag_or_versor_calibrated, voltAnchor->resonanceFrequency[2]);
            float V_rx_4 = V_from_B(B_field_vector_4, tag_or_versor_calibrated, voltAnchor->resonanceFrequency[3]);
            PredictedVoltages[0] = V_rx_1;
            PredictedVoltages[1] = V_rx_2;
            PredictedVoltages[2] = V_rx_3;
            PredictedVoltages[3] = V_rx_4;
            DEBUG_PRINT("PredictedVoltages[0] = %f\n", (double)PredictedVoltages[0]);
            DEBUG_PRINT("PredictedVoltages[1] = %f\n", (double)PredictedVoltages[1]);
            DEBUG_PRINT("PredictedVoltages[2] = %f\n", (double)PredictedVoltages[2]);
            DEBUG_PRINT("PredictedVoltages[3] = %f\n", (double)PredictedVoltages[3]);

            CG_a1 = meanData_a1 / PredictedVoltages[0];
            CG_a2 = meanData_a2 / PredictedVoltages[1];
            CG_a3 = meanData_a3 / PredictedVoltages[2];
            CG_a4 = meanData_a4 / PredictedVoltages[3];
            DEBUG_PRINT("CG_a1 = %f\n", (double)CG_a1);
            DEBUG_PRINT("CG_a2 = %f\n", (double)CG_a2);
            DEBUG_PRINT("CG_a3 = %f\n", (double)CG_a3);
            DEBUG_PRINT("CG_a4 = %f\n", (double)CG_a4);
        }
        else
        {
            // start_timer();     // start the timer.
            // it1 = get_timer(); // store current cycle-count in a local

            // computing the B field for each of the 4 anchors
            float tag_pos_predicted[3] = {this->S[KC_STATE_X], this->S[KC_STATE_Y], this->S[KC_STATE_Z]};

            T_pre_x = tag_pos_predicted[0];
            T_pre_y = tag_pos_predicted[1];
            T_pre_z = tag_pos_predicted[2];

            // sarebbe il prodotto tra la matrice di rotazione e il versore  [0,0,1] iniziale
            // float tag_or_versor[3] = {this->R[0][2], this->R[1][2], this->R[2][2]};

            float RotationMatrix[3][3];

            estimatorKalmanGetEstimatedRot((float *)RotationMatrix);
            float tag_or_versor[3] = {RotationMatrix[0][2], RotationMatrix[1][2], RotationMatrix[2][2]};

            T_orientation[0] = tag_or_versor[0];
            T_orientation[1] = tag_or_versor[1];
            T_orientation[2] = tag_or_versor[2];
            // // SEMBRA OK [0,0,1]
            // // if (pos_idx % 10 == 0)
            // // {
            // //     DEBUG_PRINT("tag_or_versor = %f,", tag_or_versor[0]);
            // //     DEBUG_PRINT("%f,", tag_or_versor[1]);
            // //     DEBUG_PRINT("%f\n", tag_or_versor[2]);
            // // }

            float anchor_1_pose[3] = {voltAnchor->x[0], voltAnchor->y[0], voltAnchor->z[0]};
            float anchor_2_pose[3] = {voltAnchor->x[1], voltAnchor->y[1], voltAnchor->z[1]};
            float anchor_3_pose[3] = {voltAnchor->x[2], voltAnchor->y[2], voltAnchor->z[2]};
            float anchor_4_pose[3] = {voltAnchor->x[3], voltAnchor->y[3], voltAnchor->z[3]};

            // ########################################### DEBUG
            // anchor_1_pose[0] = 0.5;
            // anchor_1_pose[1] = 0.5;
            // anchor_1_pose[2] = 0;

            // anchor_2_pose[0] = -0.5;
            // anchor_2_pose[1] = 0.5;
            // anchor_2_pose[2] = 0;

            // anchor_3_pose[0] = 0.5;
            // anchor_3_pose[1] = -0.5;
            // anchor_3_pose[2] = 0;

            // anchor_4_pose[0] = -0.5;
            // anchor_4_pose[1] = -0.5;
            // anchor_4_pose[2] = 0;

            // voltAnchor->x[0] = anchor_1_pose[0];
            // voltAnchor->y[0] = anchor_1_pose[1];
            // voltAnchor->z[0] = anchor_1_pose[2];

            // voltAnchor->x[1] = anchor_2_pose[0];
            // voltAnchor->y[1] = anchor_2_pose[1];
            // voltAnchor->z[1] = anchor_2_pose[2];

            // voltAnchor->x[2] = anchor_3_pose[0];
            // voltAnchor->y[2] = anchor_3_pose[1];
            // voltAnchor->z[2] = anchor_3_pose[2];

            // voltAnchor->x[3] = anchor_4_pose[0];
            // voltAnchor->y[3] = anchor_4_pose[1];
            // voltAnchor->z[3] = anchor_4_pose[2];

            // // CHECK THE ASSIGNEMENT WORKED
            // // DEBUG_PRINT("voltAnchor->x[0] = %f\n", voltAnchor->x[0]);

            // // add noise on the measurments
            // voltAnchor->measuredVolt[0] = Volts_500_4_lut[pos_idx][0]; // + generate_gaussian_noise(0.0, voltAnchor->stdDev[0]);
            // voltAnchor->measuredVolt[1] = Volts_500_4_lut[pos_idx][1]; // + generate_gaussian_noise(0.0, voltAnchor->stdDev[1]);
            // voltAnchor->measuredVolt[2] = Volts_500_4_lut[pos_idx][2]; // + generate_gaussian_noise(0.0, voltAnchor->stdDev[2]);
            // voltAnchor->measuredVolt[3] = Volts_500_4_lut[pos_idx][3]; // + generate_gaussian_noise(0.0, voltAnchor->stdDev[3]);

            // // DEBUG_PRINT("voltAnchor->measuredVolt[0] = %f\n", voltAnchor->measuredVolt[0]);
            // pos_idx = pos_idx + 1;
            // // if (pos_idx % 10 == 0)
            // // {
            // //     DEBUG_PRINT("pos_idx = %d\n", pos_idx);
            // // }
            // // ########################################### DEBUG

            MeasuredVoltages_raw[0] = voltAnchor->measuredVolt[0];
            MeasuredVoltages_raw[1] = voltAnchor->measuredVolt[1];
            MeasuredVoltages_raw[2] = voltAnchor->measuredVolt[2];
            MeasuredVoltages_raw[3] = voltAnchor->measuredVolt[3];

            MeasuredVoltages_calibrated[0] = voltAnchor->measuredVolt[0] / CG_a1;
            MeasuredVoltages_calibrated[1] = voltAnchor->measuredVolt[1] / CG_a2;
            MeasuredVoltages_calibrated[2] = voltAnchor->measuredVolt[2] / CG_a3;
            MeasuredVoltages_calibrated[3] = voltAnchor->measuredVolt[3] / CG_a4;

            get_B_field_for_a_Anchor(anchor_1_pose, tag_pos_predicted, tag_or_versor, B_field_vector_1);
            get_B_field_for_a_Anchor(anchor_2_pose, tag_pos_predicted, tag_or_versor, B_field_vector_2);
            get_B_field_for_a_Anchor(anchor_3_pose, tag_pos_predicted, tag_or_versor, B_field_vector_3);
            get_B_field_for_a_Anchor(anchor_4_pose, tag_pos_predicted, tag_or_versor, B_field_vector_4);

            B[0][0] = B_field_vector_1[0];
            B[0][1] = B_field_vector_1[1];
            B[0][2] = B_field_vector_1[2];
            // if (pos_idx % 10 == 0)
            // {
            //     DEBUG_PRINT("B_ANCHOR_0 = %f,", B[0][0] * 10e12);
            //     DEBUG_PRINT("%f,", B[0][1] * 10e12);
            //     DEBUG_PRINT("%f\n", B[0][2] * 10e10);
            // }

            B[1][0] = B_field_vector_2[0];
            B[1][1] = B_field_vector_2[1];
            B[1][2] = B_field_vector_2[2];

            B[2][0] = B_field_vector_3[0];
            B[2][1] = B_field_vector_3[1];
            B[2][2] = B_field_vector_3[2];

            B[3][0] = B_field_vector_4[0];
            B[3][1] = B_field_vector_4[1];
            B[3][2] = B_field_vector_4[2];

            // computing the V_rx for each of the 4 anchors
            float V_rx_1 = V_from_B(B_field_vector_1, tag_or_versor, voltAnchor->resonanceFrequency[0]);
            float V_rx_2 = V_from_B(B_field_vector_2, tag_or_versor, voltAnchor->resonanceFrequency[1]);
            float V_rx_3 = V_from_B(B_field_vector_3, tag_or_versor, voltAnchor->resonanceFrequency[2]);
            float V_rx_4 = V_from_B(B_field_vector_4, tag_or_versor, voltAnchor->resonanceFrequency[3]);

            PredictedVoltages[0] = V_rx_1;
            PredictedVoltages[1] = V_rx_2;
            PredictedVoltages[2] = V_rx_3;
            PredictedVoltages[3] = V_rx_4;

            // computing the V_rx_derivate for each of the 4 anchors
            float V_rx_derivative_x[4];
            float V_rx_derivative_y[4];
            float V_rx_derivative_z[4];

            V_rx_derivate_x(
                tag_pos_predicted[0], tag_pos_predicted[1], tag_pos_predicted[2],
                tag_or_versor[0], tag_or_versor[1], tag_or_versor[2],
                N_WOUNDS, CURRENT, COIL_SURFACE, G_INA,
                voltAnchor->x[0], voltAnchor->y[0], voltAnchor->z[0],
                voltAnchor->x[1], voltAnchor->y[1], voltAnchor->z[1],
                voltAnchor->x[2], voltAnchor->y[2], voltAnchor->z[2],
                voltAnchor->x[3], voltAnchor->y[3], voltAnchor->z[3],
                voltAnchor->resonanceFrequency[0], voltAnchor->resonanceFrequency[1], voltAnchor->resonanceFrequency[2], voltAnchor->resonanceFrequency[3],
                V_rx_derivative_x);

            V_rx_derivate_y(
                tag_pos_predicted[0], tag_pos_predicted[1], tag_pos_predicted[2],
                tag_or_versor[0], tag_or_versor[1], tag_or_versor[2],
                N_WOUNDS, CURRENT, COIL_SURFACE, G_INA,
                voltAnchor->x[0], voltAnchor->y[0], voltAnchor->z[0],
                voltAnchor->x[1], voltAnchor->y[1], voltAnchor->z[1],
                voltAnchor->x[2], voltAnchor->y[2], voltAnchor->z[2],
                voltAnchor->x[3], voltAnchor->y[3], voltAnchor->z[3],
                voltAnchor->resonanceFrequency[0], voltAnchor->resonanceFrequency[1], voltAnchor->resonanceFrequency[2], voltAnchor->resonanceFrequency[3],
                V_rx_derivative_y);

            V_rx_derivate_z(
                tag_pos_predicted[0], tag_pos_predicted[1], tag_pos_predicted[2],
                tag_or_versor[0], tag_or_versor[1], tag_or_versor[2],
                N_WOUNDS, CURRENT, COIL_SURFACE, G_INA,
                voltAnchor->x[0], voltAnchor->y[0], voltAnchor->z[0],
                voltAnchor->x[1], voltAnchor->y[1], voltAnchor->z[1],
                voltAnchor->x[2], voltAnchor->y[2], voltAnchor->z[2],
                voltAnchor->x[3], voltAnchor->y[3], voltAnchor->z[3],
                voltAnchor->resonanceFrequency[0], voltAnchor->resonanceFrequency[1], voltAnchor->resonanceFrequency[2], voltAnchor->resonanceFrequency[3],
                V_rx_derivative_z);

            // it2 = get_timer() - it1; // Derive the cycle-count difference
            // stop_timer();

            // Debugging stuff
            for (int i = 0; i < 4; i++)
            {
                Derivative_x[i] = V_rx_derivative_x[i];
                Derivative_y[i] = V_rx_derivative_y[i];
                Derivative_z[i] = V_rx_derivative_z[i];
            }

            // create the matrix to update the kalman matrix
            float h_1[KC_STATE_DIM] = {0};
            arm_matrix_instance_f32 H_1 = {1, KC_STATE_DIM, h_1};
            h_1[KC_STATE_X] = V_rx_derivative_x[1];
            h_1[KC_STATE_Y] = V_rx_derivative_y[1];
            // h_1[KC_STATE_Z] = V_rx_derivative_z[1];

            float h_2[KC_STATE_DIM] = {0};
            arm_matrix_instance_f32 H_2 = {1, KC_STATE_DIM, h_2};
            h_2[KC_STATE_X] = V_rx_derivative_x[2];
            h_2[KC_STATE_Y] = V_rx_derivative_y[2];
            // h_2[KC_STATE_Z] = V_rx_derivative_z[2];

            float h_3[KC_STATE_DIM] = {0};
            arm_matrix_instance_f32 H_3 = {1, KC_STATE_DIM, h_3};
            h_3[KC_STATE_X] = V_rx_derivative_x[3];
            h_3[KC_STATE_Y] = V_rx_derivative_y[3];
            h_3[KC_STATE_Z] = V_rx_derivative_z[3];

            float h_4[KC_STATE_DIM] = {0};
            arm_matrix_instance_f32 H_4 = {1, KC_STATE_DIM, h_4};
            h_4[KC_STATE_X] = V_rx_derivative_x[4];
            h_4[KC_STATE_Y] = V_rx_derivative_y[4];
            h_4[KC_STATE_Z] = V_rx_derivative_z[4];

            float error_anchor1 = MeasuredVoltages_calibrated[0] - V_rx_1;
            float error_anchor2 = MeasuredVoltages_calibrated[1] - V_rx_2;
            float error_anchor3 = MeasuredVoltages_calibrated[2] - V_rx_3;
            float error_anchor4 = MeasuredVoltages_calibrated[3] - V_rx_4;
            E_1 = error_anchor1;
            E_2 = error_anchor2;
            E_3 = error_anchor3;
            E_4 = error_anchor4;

            CalibrationRapport_anchor1 = MeasuredVoltages_calibrated[0] / PredictedVoltages[0];
            CalibrationRapport_anchor2 = MeasuredVoltages_calibrated[1] / PredictedVoltages[1];
            CalibrationRapport_anchor3 = MeasuredVoltages_calibrated[2] / PredictedVoltages[2];
            CalibrationRapport_anchor4 = MeasuredVoltages_calibrated[3] / PredictedVoltages[3];

            // Potrebbe essere troppo piccolo l'errore. o le stime troppo diverse?!?!?!
            // kalmanCoreScalarUpdate(this, &H_1, error_anchor1, voltAnchor->stdDev[0]);
            // kalmanCoreScalarUpdate(this, &H_2, error_anchor2, voltAnchor->stdDev[1]);
            // kalmanCoreScalarUpdate(this, &H_3, error_anchor3, voltAnchor->stdDev[2]);
            // kalmanCoreScalarUpdate(this, &H_4, error_anchor4, voltAnchor->stdDev[3]);
        }
    }
}

LOG_GROUP_START(Dipole_Model)

// LOG_ADD(LOG_UINT32, CPUCycle, &it2)

LOG_ADD(LOG_FLOAT, D_x_0, &Derivative_x[0])
LOG_ADD(LOG_FLOAT, D_x_1, &Derivative_x[1])
LOG_ADD(LOG_FLOAT, D_x_2, &Derivative_x[2])
LOG_ADD(LOG_FLOAT, D_x_3, &Derivative_x[3])
LOG_ADD(LOG_FLOAT, D_y_0, &Derivative_y[0])
LOG_ADD(LOG_FLOAT, D_y_1, &Derivative_y[1])
LOG_ADD(LOG_FLOAT, D_y_2, &Derivative_y[2])
LOG_ADD(LOG_FLOAT, D_y_3, &Derivative_y[3])
LOG_ADD(LOG_FLOAT, D_z_0, &Derivative_z[0])
LOG_ADD(LOG_FLOAT, D_z_1, &Derivative_z[1])
LOG_ADD(LOG_FLOAT, D_z_2, &Derivative_z[2])
LOG_ADD(LOG_FLOAT, D_z_3, &Derivative_z[3])

LOG_ADD(LOG_FLOAT, CG_A1, &CG_a1)
LOG_ADD(LOG_FLOAT, CG_A2, &CG_a2)
LOG_ADD(LOG_FLOAT, CG_A3, &CG_a3)
LOG_ADD(LOG_FLOAT, CG_A4, &CG_a4)

LOG_ADD(LOG_FLOAT, C_R_1, &CalibrationRapport_anchor1)
LOG_ADD(LOG_FLOAT, C_R_2, &CalibrationRapport_anchor2)
LOG_ADD(LOG_FLOAT, C_R_3, &CalibrationRapport_anchor3)
LOG_ADD(LOG_FLOAT, C_R_4, &CalibrationRapport_anchor4)

LOG_ADD(LOG_FLOAT, T_pred_x, &T_pre_x)
LOG_ADD(LOG_FLOAT, T_pred_y, &T_pre_y)
LOG_ADD(LOG_FLOAT, T_pred_z, &T_pre_z)

LOG_ADD(LOG_FLOAT, T_or_0, &T_orientation[0])
LOG_ADD(LOG_FLOAT, T_or_1, &T_orientation[1])
LOG_ADD(LOG_FLOAT, T_or_2, &T_orientation[2])

LOG_ADD(LOG_FLOAT, Er_1, &E_1)
LOG_ADD(LOG_FLOAT, Er_2, &E_2)
LOG_ADD(LOG_FLOAT, Er_3, &E_3)
LOG_ADD(LOG_FLOAT, Er_4, &E_4)

LOG_ADD(LOG_FLOAT, MC_V_0, &MeasuredVoltages_calibrated[0])
LOG_ADD(LOG_FLOAT, MC_V_1, &MeasuredVoltages_calibrated[1])
LOG_ADD(LOG_FLOAT, MC_V_2, &MeasuredVoltages_calibrated[2])
LOG_ADD(LOG_FLOAT, MC_V_3, &MeasuredVoltages_calibrated[3])

LOG_ADD(LOG_FLOAT, MR_V_0, &MeasuredVoltages_raw[0])
LOG_ADD(LOG_FLOAT, MR_V_1, &MeasuredVoltages_raw[1])
LOG_ADD(LOG_FLOAT, MR_V_2, &MeasuredVoltages_raw[2])
LOG_ADD(LOG_FLOAT, MR_V_3, &MeasuredVoltages_raw[3])

LOG_ADD(LOG_FLOAT, P_V_0, &PredictedVoltages[0])
LOG_ADD(LOG_FLOAT, P_V_1, &PredictedVoltages[1])
LOG_ADD(LOG_FLOAT, P_V_2, &PredictedVoltages[2])
LOG_ADD(LOG_FLOAT, P_V_3, &PredictedVoltages[3])

LOG_ADD(LOG_FLOAT, B_0_0, &B[0][0])
LOG_ADD(LOG_FLOAT, B_0_1, &B[0][1])
LOG_ADD(LOG_FLOAT, B_0_2, &B[0][2])

LOG_ADD(LOG_FLOAT, B_1_0, &B[1][0])
LOG_ADD(LOG_FLOAT, B_1_1, &B[1][1])
LOG_ADD(LOG_FLOAT, B_1_2, &B[1][2])

LOG_ADD(LOG_FLOAT, B_2_0, &B[2][0])
LOG_ADD(LOG_FLOAT, B_2_1, &B[2][1])
LOG_ADD(LOG_FLOAT, B_2_2, &B[2][2])

LOG_ADD(LOG_FLOAT, B_3_0, &B[3][0])
LOG_ADD(LOG_FLOAT, B_3_1, &B[3][1])
LOG_ADD(LOG_FLOAT, B_3_2, &B[3][2])

LOG_GROUP_STOP(Dipole_Model)

PARAM_GROUP_START(Dipole_Params)
PARAM_ADD(PARAM_UINT16, calibTic, &currentCalibrationTick)
PARAM_GROUP_STOP(Dipole_Params)