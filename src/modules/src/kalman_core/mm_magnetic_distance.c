#include "math.h"
#include "mm_magnetic_distance.h"
#include "log.h"

#define G_INA 1000.0f
#define RAY 0.019f
#define N_WOUNDS 5.0f
#define COIL_SURFACE (RAY * RAY * PI)
#define CURRENT 1.5f
#define MU_0 1.25663706e-6f

#define start_timer() *((volatile uint32_t *)0xE0001000) = 0x40000001 // Enable CYCCNT register
#define stop_timer() *((volatile uint32_t *)0xE0001000) = 0x40000000  // Disable CYCCNT register
#define get_timer() *((volatile uint32_t *)0xE0001004)

uint32_t it1, it2; // start and stop flag
float32_t MeasuredVoltages[4] = {0.0f, 0.0f, 0.0f, 0.0f};
float32_t PredictedVoltages[4] = {0.0f, 0.0f, 0.0f, 0.0f};
float32_t Derivative_x[4] = {0.0f, 0.0f, 0.0f, 0.0f};
float32_t Derivative_y[4] = {0.0f, 0.0f, 0.0f, 0.0f};
float32_t Derivative_z[4] = {0.0f, 0.0f, 0.0f, 0.0f};
float32_t B[4][3] = {{0.0f, 0.0f, 0.0f},
                     {0.0f, 0.0f, 0.0f},
                     {0.0f, 0.0f, 0.0f},
                     {0.0f, 0.0f, 0.0f}};

// float sign(float x)
// {
//     if (x > 0)
//         return 1.0f;
//     else if (x < 0)
//         return -1.0f;
//     else
//         return 0.0f;
// }

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
    float sum = 0.0;
    for (int i = 0; i < length; i++)
    {
        float diff = a[i] - b[i];
        sum += diff * diff;
    }
    return (float)sqrtf(sum);
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
                           copysign(1.0, MU_0 * W_x * t7 * t51 * (t12 - t98) * (-1.0f / 4.0f) -
                                             (MU_0 * W_y * t7 * t51 * (t13 - t102)) / 4.0f -
                                             (MU_0 * W_z * t7 * t51 * (t14 - t106)) / 4.0f) *
                           ((MU_0 * W_y * t7 * t51 * (t19 * t50 * t122 + (t19 * t30 * t51 * t94) / 2.0f)) / 4.0f +
                            (MU_0 * W_z * t7 * t51 * (t23 * t50 * t122 + (t23 * t30 * t51 * t94) / 2.0f)) / 4.0f +
                            (MU_0 * W_x * t7 * t51 * (-t50 * t94 + t15 * t50 * t122 + (t15 * t30 * t51 * t94) / 2.0f)) / 4.0f -
                            MU_0 * W_x * t7 * t30 * t53 * (t12 - t98) * (3.0f / 8.0f) -
                            MU_0 * W_y * t7 * t30 * t53 * (t13 - t102) * (3.0f / 8.0f) -
                            MU_0 * W_z * t7 * t30 * t53 * (t14 - t106) * (3.0f / 8.0f)) *
                           2.0;

    float t107 = t24 * t52 * t95;
    V_rx_derivative_x[1] = F_1 * G_ina * N * S * PI *
                           copysign(1.0, MU_0 * W_x * t7 * t54 * (t12 - t99) * (-1.0f / 4.0f) -
                                             (MU_0 * W_y * t7 * t54 * (t13 - t103)) / 4.0f -
                                             (MU_0 * W_z * t7 * t54 * (t14 - t107)) / 4.0f) *
                           ((MU_0 * W_y * t7 * t54 * (t20 * t52 * t123 + (t20 * t31 * t54 * t95) / 2.0f)) / 4.0f +
                            (MU_0 * W_z * t7 * t54 * (t24 * t52 * t123 + (t24 * t31 * t54 * t95) / 2.0f)) / 4.0f +
                            (MU_0 * W_x * t7 * t54 * (-t52 * t95 + t16 * t52 * t123 + (t16 * t31 * t54 * t95) / 2.0f)) / 4.0f -
                            MU_0 * W_x * t7 * t31 * t56 * (t12 - t99) * (3.0f / 8.0f) -
                            MU_0 * W_y * t7 * t31 * t56 * (t13 - t103) * (3.0f / 8.0f) -
                            MU_0 * W_z * t7 * t31 * t56 * (t14 - t107) * (3.0f / 8.0f)) *
                           2.0;

    float t108 = t25 * t55 * t96;
    V_rx_derivative_x[2] = F_2 * G_ina * N * S * PI *
                           copysign(1.0, MU_0 * W_x * t7 * t57 * (t12 - t100) * (-1.0f / 4.0f) -
                                             (MU_0 * W_y * t7 * t57 * (t13 - t104)) / 4.0f -
                                             (MU_0 * W_z * t7 * t57 * (t14 - t108)) / 4.0f) *
                           ((MU_0 * W_y * t7 * t57 * (t21 * t55 * t124 + (t21 * t32 * t57 * t96) / 2.0f)) / 4.0f +
                            (MU_0 * W_z * t7 * t57 * (t25 * t55 * t124 + (t25 * t32 * t57 * t96) / 2.0f)) / 4.0f +
                            (MU_0 * W_x * t7 * t57 * (-t55 * t96 + t17 * t55 * t124 + (t17 * t32 * t57 * t96) / 2.0f)) / 4.0f -
                            MU_0 * W_x * t7 * t32 * t59 * (t12 - t100) * (3.0f / 8.0f) -
                            MU_0 * W_y * t7 * t32 * t59 * (t13 - t104) * (3.0f / 8.0f) -
                            MU_0 * W_z * t7 * t32 * t59 * (t14 - t108) * (3.0f / 8.0f)) *
                           2.0;

    float t109 = t26 * t58 * t97;
    V_rx_derivative_x[3] = F_3 * G_ina * N * S * PI *
                           copysign(1.0, MU_0 * W_x * t7 * t60 * (t12 - t101) * (-1.0f / 4.0f) -
                                             (MU_0 * W_y * t7 * t60 * (t13 - t105)) / 4.0f -
                                             (MU_0 * W_z * t7 * t60 * (t14 - t109)) / 4.0f) *
                           ((MU_0 * W_y * t7 * t60 * (t22 * t58 * t125 + (t22 * t33 * t60 * t97) / 2.0f)) / 4.0f +
                            (MU_0 * W_z * t7 * t60 * (t26 * t58 * t125 + (t26 * t33 * t60 * t97) / 2.0f)) / 4.0f +
                            (MU_0 * W_x * t7 * t60 * (-t58 * t97 + t18 * t58 * t125 + (t18 * t33 * t60 * t97) / 2.0f)) / 4.0f -
                            MU_0 * W_x * t7 * t33 * t61 * (t12 - t101) * (3.0f / 8.0f) -
                            MU_0 * W_y * t7 * t33 * t61 * (t13 - t105) * (3.0f / 8.0f) -
                            MU_0 * W_z * t7 * t33 * t61 * (t14 - t109) * (3.0f / 8.0f)) *
                           2.0;
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

float get_B_field_for_a_Anchor(float *anchor_pos,
                               float *tag_pos,
                               float *tag_or_versor,
                               float *B_field)
{
    float tx_rx_versor[3];
    getversor(anchor_pos, tag_pos, tx_rx_versor, 3);

    float magnetic_dipole_moment_tx[3] = {
        N_WOUNDS * COIL_SURFACE * CURRENT * tag_or_versor[0],
        N_WOUNDS * COIL_SURFACE * CURRENT * tag_or_versor[1],
        N_WOUNDS * COIL_SURFACE * CURRENT * tag_or_versor[2]};
    float tx_rx_distance = euclidean_distance(anchor_pos, tag_pos, 3);

    float dot_product_B_temp = dot_product(magnetic_dipole_moment_tx, tx_rx_versor, 3);

    float constant_Bfield_constant_1 = (MU_0 / (4 * PI)) / powf(tx_rx_distance, 3);
    float B_temp[3] = {
        constant_Bfield_constant_1 * (3 * dot_product_B_temp * tx_rx_versor[0] - magnetic_dipole_moment_tx[0]),
        constant_Bfield_constant_1 * (3 * dot_product_B_temp * tx_rx_versor[1] - magnetic_dipole_moment_tx[1]),
        constant_Bfield_constant_1 * (3 * dot_product_B_temp * tx_rx_versor[2] - magnetic_dipole_moment_tx[2])};

    float magnetic_dipole_moment_tx_magnitude = euclidean_distance(magnetic_dipole_moment_tx, magnetic_dipole_moment_tx, 3);

    B_field[0] = B_temp[0] * magnetic_dipole_moment_tx_magnitude,
    B_field[1] = B_temp[1] * magnetic_dipole_moment_tx_magnitude,
    B_field[2] = B_temp[2] * magnetic_dipole_moment_tx_magnitude;
}

float V_from_B(float *B_field, float *rx_versor, int resonanceFreq)
{

    float dot_product_V = dot_product(B_field, rx_versor, 3);
    float V = G_INA * fabs(2 * M_PI * resonanceFreq * COIL_SURFACE * N_WOUNDS * dot_product_V);
    return V;
}

float kalmanCoreUpdateWithVolt(kalmanCoreData_t *this,
                               voltMeasurement_t *voltAnchor)
{

    start_timer();     // start the timer.
    it1 = get_timer(); // store current cycle-count in a local
    MeasuredVoltages[0] = voltAnchor->measuredVolt[0];
    MeasuredVoltages[1] = voltAnchor->measuredVolt[1];
    MeasuredVoltages[2] = voltAnchor->measuredVolt[2];
    MeasuredVoltages[3] = voltAnchor->measuredVolt[3];

    // computing the B field for each of the 4 anchors
    float tag_pos_predicted[3] = {this->S[KC_STATE_X], this->S[KC_STATE_Y], this->S[KC_STATE_Z]};

    // sarebbe il prodotto tra la matrice di rotazione e il versore  [0,0,1] iniziale
    float tag_or_versor[3] = {this->R[0][2], this->R[1][2], this->R[2][2]};

    float B_field_vector_1[3];
    float B_field_vector_2[3];
    float B_field_vector_3[3];
    float B_field_vector_4[3];

    float anchor_1_pose[3] = {voltAnchor->x[0], voltAnchor->y[0], voltAnchor->z[0]};
    float anchor_2_pose[3] = {voltAnchor->x[1], voltAnchor->y[1], voltAnchor->z[1]};
    float anchor_3_pose[3] = {voltAnchor->x[2], voltAnchor->y[2], voltAnchor->z[2]};
    float anchor_4_pose[3] = {voltAnchor->x[3], voltAnchor->y[3], voltAnchor->z[3]};

    get_B_field_for_a_Anchor(anchor_1_pose, tag_pos_predicted, tag_or_versor, B_field_vector_1);
    get_B_field_for_a_Anchor(anchor_2_pose, tag_pos_predicted, tag_or_versor, B_field_vector_2);
    get_B_field_for_a_Anchor(anchor_3_pose, tag_pos_predicted, tag_or_versor, B_field_vector_3);
    get_B_field_for_a_Anchor(anchor_4_pose, tag_pos_predicted, tag_or_versor, B_field_vector_4);

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

    it2 = get_timer() - it1; // Derive the cycle-count difference
    stop_timer();

    // create the matrix to update the kalman matrix
    float h_1[KC_STATE_DIM] = {0};
    arm_matrix_instance_f32 H_1 = {1, KC_STATE_DIM, h_1};
    h_1[KC_STATE_X] = V_rx_derivative_x[1];
    h_1[KC_STATE_Y] = V_rx_derivative_y[1];
    h_1[KC_STATE_Z] = V_rx_derivative_z[1];

    float h_2[KC_STATE_DIM] = {0};
    arm_matrix_instance_f32 H_2 = {1, KC_STATE_DIM, h_2};
    h_2[KC_STATE_X] = V_rx_derivative_x[2];
    h_2[KC_STATE_Y] = V_rx_derivative_y[2];
    h_2[KC_STATE_Z] = V_rx_derivative_z[2];

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
    // Potrebbe essere troppo piccolo l'errore. o le stime troppo diverse?!?!?!
    kalmanCoreScalarUpdate(this, &H_1, voltAnchor->measuredVolt[0] - V_rx_1, voltAnchor->stdDev[0]);
    kalmanCoreScalarUpdate(this, &H_2, voltAnchor->measuredVolt[1] - V_rx_2, voltAnchor->stdDev[1]);
    kalmanCoreScalarUpdate(this, &H_3, voltAnchor->measuredVolt[2] - V_rx_3, voltAnchor->stdDev[2]);
    kalmanCoreScalarUpdate(this, &H_4, voltAnchor->measuredVolt[3] - V_rx_4, voltAnchor->stdDev[3]);
}

LOG_GROUP_START(Dipole_Model)

// LOG_ADD(LOG_UINT32, CPUCycle, &it2)

LOG_ADD(LOG_FLOAT, M_V_0, &MeasuredVoltages[0])
LOG_ADD(LOG_FLOAT, M_V_1, &MeasuredVoltages[1])
LOG_ADD(LOG_FLOAT, M_V_2, &MeasuredVoltages[2])
LOG_ADD(LOG_FLOAT, M_V_3, &MeasuredVoltages[3])

LOG_ADD(LOG_FLOAT, P_V_0, &PredictedVoltages[0])
LOG_ADD(LOG_FLOAT, P_V_1, &PredictedVoltages[1])
LOG_ADD(LOG_FLOAT, P_V_2, &PredictedVoltages[2])
LOG_ADD(LOG_FLOAT, P_V_3, &PredictedVoltages[3])

LOG_ADD(LOG_FLOAT, B_0_0, &B[0][0])
LOG_ADD(LOG_FLOAT, B_0_1, &B[0][1])
LOG_ADD(LOG_FLOAT, B_0_2, &B[0][2])
LOG_ADD(LOG_FLOAT, B_0_3, &B[0][3])

LOG_ADD(LOG_FLOAT, B_1_0, &B[1][0])
LOG_ADD(LOG_FLOAT, B_1_1, &B[1][1])
LOG_ADD(LOG_FLOAT, B_1_2, &B[1][2])
LOG_ADD(LOG_FLOAT, B_1_3, &B[1][3])

LOG_ADD(LOG_FLOAT, B_2_0, &B[2][0])
LOG_ADD(LOG_FLOAT, B_2_1, &B[2][1])
LOG_ADD(LOG_FLOAT, B_2_2, &B[2][2])
LOG_ADD(LOG_FLOAT, B_2_3, &B[2][3])

LOG_ADD(LOG_FLOAT, B_3_0, &B[3][0])
LOG_ADD(LOG_FLOAT, B_3_1, &B[3][1])
LOG_ADD(LOG_FLOAT, B_3_2, &B[3][2])
LOG_ADD(LOG_FLOAT, B_3_3, &B[3][3])

LOG_GROUP_STOP(Dipole_Model)
