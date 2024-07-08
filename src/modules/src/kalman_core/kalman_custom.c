#include "kalman_custom.h"
#include <string.h>
#include <stdint.h>
#include "debug.h"

void matmul_3x3x3(float A[3][3], float B[3][3], float R[3][3])
{
    float A_copy[3][3];
    float B_copy[3][3];
    memcpy(A_copy, A, 9 * sizeof(float));
    memcpy(B_copy, B, 9 * sizeof(float));

    for (int16_t i = 0; i < 3; i++)
    {
        for (int16_t j = 0; j < 3; j++)
        {
            R[i][j] = 0;
            for (int16_t k = 0; k < 3; k++)
                R[i][j] += A_copy[i][k] * B_copy[k][j];
        }
    }
}

void transpose_matrix_3x3(float original[3][3], float transpose[3][3])
{
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            transpose[j][i] = original[i][j];
        }
    }
}

void matmul_1x3x3(float vector[3], float matrix[3][3], float result[3])
{
    for (int i = 0; i < 3; i++)
    {
        result[i] = 0; // Inizializza l'elemento del risultato
        for (int j = 0; j < 3; j++)
        {
            result[i] += vector[j] * matrix[j][i]; // Moltiplicazione e somma
        }
    }
}

void matinv_3x3(float A[3][3], float A_inv[3][3])
{
    float det = 0.0f;
    int i, j;
    for (i = 0; i < 3; i++)
        det = det + (A[0][i] * (A[1][(i + 1) % 3] * A[2][(i + 2) % 3] -
                                A[1][(i + 2) % 3] * A[2][(i + 1) % 3]));

    float A_copy[3][3];
    memcpy(A_copy, A, 9 * sizeof(float));

    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            A_inv[j][i] = ((A_copy[(i + 1) % 3][(j + 1) % 3] *
                            A_copy[(i + 2) % 3][(j + 2) % 3]) -
                           (A_copy[(i + 1) % 3][(j + 2) % 3] *
                            A_copy[(i + 2) % 3][(j + 1) % 3])) /
                          det;
}

void matmul_3x3x1(float matrix[3][3], float array[3], float result[3])
{
    for (int i = 0; i < 3; i++)
    {
        result[i] = 0; // Inizializza l'elemento del risultato
        for (int j = 0; j < 3; j++)
        {
            result[i] += matrix[i][j] * array[j]; // Moltiplicazione e somma
        }
    }
}

float dot_prod(float *a, float *b)
{
    float result = 0;
    for (int i = 0; i < 3; i++)
    {
        result += a[i] * b[i];
    }
    return result;
}

void multiply_by_scalar(float *a, float scalar)
{
    for (int i = 0; i < 3; i++)
    {
        a[i] *= scalar;
    }
}

void multiply_3x1_with_1x3(float vector3x1[3], float vector1x3[3], float resultMatrix[3][3])
{
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            resultMatrix[i][j] = vector3x1[i] * vector1x3[j];
        }
    }
}

void multiply_matrix_3x3_by_scalar(float matrix[3][3], float scalar)
{
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            matrix[i][j] = matrix[i][j] * scalar;
        }
    }
}

void kalman_custom_predict(ekf_t *ekf)
{
    for (int i = 0; i < 3; i++)
    {
        ekf->P[i][i] += ekf->Q;
    }
}

void kalman_custom_update(ekf_t *ekf, float H[3], float delta)
{
    // Compute Kalman gain K
    float K[3];
    matmul_1x3x3(H, ekf->P, K);
    float inverse = dot_prod(K, H);
    inverse += ekf->R;

    DEBUG_PRINT("inverse: %f\n", inverse);
    DEBUG_PRINT("P: %f %f %f\n", ekf->P[0][0], ekf->P[1][1], ekf->P[2][2]);

    if (fabs(inverse) > 0.00001f)
    {
        inverse = 1 / inverse;
    }
    else
    {
        DEBUG_PRINT("Inverse is 0\n");
        return;
    }
    matmul_3x3x1(ekf->P, H, K);
    multiply_by_scalar(K, inverse);

    DEBUG_PRINT("K: %f %f %f\n", K[0], K[1], K[2]);

    // Compute x
    float correction[3]; // correction = K * delta
    correction[0] = K[0] * delta;
    correction[1] = K[1] * delta;
    correction[2] = K[2] * delta;

    ekf->x[0] += correction[0];
    ekf->x[1] += correction[1];
    ekf->x[2] += correction[2];

    // Compute P
    float IKH[3][3];
    multiply_3x1_with_1x3(K, H, IKH);
    IKH[0][0] = 1 - IKH[0][0];
    IKH[1][1] = 1 - IKH[1][1];
    IKH[2][2] = 1 - IKH[2][2];

    float IKH_T[3][3];
    transpose_matrix_3x3(IKH, IKH_T);

    float Pm[3][3];
    matmul_3x3x3(IKH, ekf->P, Pm);
    matmul_3x3x3(Pm, IKH_T, Pm);

    float KRKT[3][3];
    multiply_3x1_with_1x3(K, K, KRKT);
    multiply_matrix_3x3_by_scalar(KRKT, ekf->R);

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            ekf->P[i][j] = Pm[i][j] + KRKT[i][j];
        }
    }
}