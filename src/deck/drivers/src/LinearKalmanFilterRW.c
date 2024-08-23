#include "usec_time.h"
#include <time.h>
#include "arm_math.h" // Libreria CMSIS-DSP per operazioni con matrici
#include "MagneticDeck.h"
#include "debug.h"

float diff_in_ms = 0;
void add_regularization(float (*matrix)[MEASURE_DIM], int dim)
{
    for (int i = 0; i < dim; i++)
    {
        matrix[i][i] += EPSILON;
    }
}

// Moltiplicazione di due matrici
void matrix_multiply(int n, int m, int p, float A[n][m], float B[m][p], float result[n][p])
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < p; j++)
        {
            result[i][j] = 0;
            for (int k = 0; k < m; k++)
            {
                result[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

// Somma di due matrici
void matrix_add(int n, int m, float A[n][m], float B[n][m], float result[n][m])
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            result[i][j] = A[i][j] + B[i][j];
        }
    }
}

// Sottrazione di due matrici
void matrix_subtract(int n, int m, float A[n][m], float B[n][m], float result[n][m])
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            result[i][j] = A[i][j] - B[i][j];
        }
    }
}

int matrix_inversion(int n, float A[n][n], float inverse[n][n])
{
    int i, j, k;
    float temp;
    float augmented[n][2 * n];

    // Creazione della matrice estesa [A | I]
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            augmented[i][j] = A[i][j];
            augmented[i][j + n] = (i == j) ? 1 : 0;
        }
    }

    // Gauss-Jordan
    for (i = 0; i < n; i++)
    {
        temp = augmented[i][i];
        if (temp == 0)

            return 0; // Matrice non invertibile
        for (j = 0; j < 2 * n; j++)
            augmented[i][j] /= temp;

        for (j = 0; j < n; j++)
        {
            if (j != i)
            {
                temp = augmented[j][i];
                for (k = 0; k < 2 * n; k++)
                    augmented[j][k] -= augmented[i][k] * temp;
            }
        }
    }

    // Estrazione della matrice inversa
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            inverse[i][j] = augmented[i][j + n];
        }
    }

    return 1; // Successo
}

// #define STATE_DIM 3   // Stato: [x, y, z]
// #define MEASURE_DIM 3 // Misure: [x, y, z]

void kalman_init(KalmanFilter *kf, float *initial_guess)
{
    // Inizializza la matrice di transizione dello stato F
    float F_init[STATE_DIM][STATE_DIM] = {
        {1, 0, 0},
        {0, 1, 0},
        {0, 0, 1}};
    for (int i = 0; i < STATE_DIM; i++)
    {
        for (int j = 0; j < STATE_DIM; j++)
        {
            kf->F[i][j] = F_init[i][j];
        }
    }

    // Inizializza la matrice di osservazione H
    float H_init[MEASURE_DIM][STATE_DIM] = {
        {1, 0, 0},
        {0, 1, 0},
        {0, 0, 1}};
    for (int i = 0; i < MEASURE_DIM; i++)
    {
        for (int j = 0; j < STATE_DIM; j++)
        {
            kf->H[i][j] = H_init[i][j];
        }
    }

    // Inizializza la matrice di covarianza del rumore di processo Q

    float Q_init[STATE_DIM][STATE_DIM] = {
        {q_kf_default * q_kf_default, 0, 0},
        {0, q_kf_default * q_kf_default, 0},
        {0, 0, q_kf_default * q_kf_default},
    };
    for (int i = 0; i < STATE_DIM; i++)
    {
        for (int j = 0; j < STATE_DIM; j++)
        {
            kf->Q[i][j] = Q_init[i][j];
        }
    }

    // Inizializza la matrice di covarianza del rumore di misura R
    float R_init[MEASURE_DIM][MEASURE_DIM] = {
        {sigma_x_kf_default * sigma_x_kf_default, 0, 0},
        {0, sigma_x_kf_default * sigma_x_kf_default, 0},
        {0, 0, sigma_x_kf_default * sigma_x_kf_default}};
    for (int i = 0; i < MEASURE_DIM; i++)
    {
        for (int j = 0; j < MEASURE_DIM; j++)
        {
            kf->R[i][j] = R_init[i][j];
        }
    }

    // Inizializza la matrice di covarianza dell'errore P
    float P_init[STATE_DIM][STATE_DIM] = {
        {1, 0, 0},
        {0, 1, 0},
        {0, 0, 1},
    };

    for (int i = 0; i < STATE_DIM; i++)
    {
        for (int j = 0; j < STATE_DIM; j++)
        {
            kf->P[i][j] = P_init[i][j];
        }
    }

    // Inizializza lo stato iniziale x_state
    for (int i = 0; i < STATE_DIM; i++)
    {
        kf->x_state[i] = initial_guess[i];
    }

    // Inizializza il vettore di misura z_meas
    for (int i = 0; i < MEASURE_DIM; i++)
    {
        kf->z_meas[i] = initial_guess[i];
    }
}

// Predizione
void kalman_predict(KalmanFilter *kf)
{
    // printf("initial guess: \n");
    // printf("x = %f, y = %f, z = %f\n", kf->x_state[0], kf->x_state[1], kf->x_state[2]);

    float temp1[STATE_DIM][STATE_DIM];
    float FT[STATE_DIM][STATE_DIM];

    // x = F * x
    for (int i = 0; i < STATE_DIM; i++)
    {
        float temp = 0;
        for (int j = 0; j < STATE_DIM; j++)
        {
            temp += kf->F[i][j] * kf->x_state[j];
        }
        kf->x_state[i] = temp;
    }

    // P = F * P * F' + Q
    // matrix_multiply(F, P, temp1);
    matrix_multiply(STATE_DIM, STATE_DIM, STATE_DIM, kf->F, kf->P, temp1);
    for (int i = 0; i < STATE_DIM; i++)
    {
        for (int j = 0; j < STATE_DIM; j++)
        {
            FT[i][j] = kf->F[j][i]; // Transpose of F
        }
    }
    // matrix_multiply(temp1, FT, P);
    matrix_multiply(STATE_DIM, STATE_DIM, STATE_DIM, temp1, FT, kf->P);

    // matrix_add(P, Q, P);
    matrix_add(STATE_DIM, STATE_DIM, kf->P, kf->Q, kf->P);
}

// Aggiornamento
void kalman_update(KalmanFilter *kf)
{
    float K[STATE_DIM][MEASURE_DIM];
    float temp1[MEASURE_DIM][STATE_DIM];
    float temp2[MEASURE_DIM][MEASURE_DIM];
    float y[MEASURE_DIM];

    float temp3[MEASURE_DIM][MEASURE_DIM];

    // y = z - H * x
    for (int i = 0; i < MEASURE_DIM; i++)
    {
        y[i] = kf->z_meas[i];
        for (int j = 0; j < STATE_DIM; j++)
        {
            y[i] -= kf->H[i][j] * kf->x_state[j];
        }
    }

    // K = P * H' * (H * P * H' + R)^-1
    // Qui puoi usare una funzione di inversione di matrice adatta alla dimensione della matrice R
    // matrix_multiply(P, H, temp1);
    // matrix_multiply(temp1, H, temp2);
    // matrix_add(temp2, R, temp2);
    matrix_multiply(STATE_DIM, MEASURE_DIM, STATE_DIM, kf->P, kf->H, temp1);
    matrix_multiply(MEASURE_DIM, STATE_DIM, MEASURE_DIM, temp1, kf->H, temp2);
    matrix_add(MEASURE_DIM, MEASURE_DIM, temp2, kf->R, temp2);

    // Inverti temp2
    // matrix_inversion(temp2, temp2);

    if (matrix_inversion(MEASURE_DIM, temp2, temp3) == 0)
    {
        DEBUG_PRINT("Matrix inversion failed\n");
        add_regularization(temp2, MEASURE_DIM);
        matrix_inversion(MEASURE_DIM, temp2, temp3);
    }

    // matrix_multiply(temp1, temp2, K);
    matrix_multiply(STATE_DIM, MEASURE_DIM, MEASURE_DIM, temp1, temp3, K);

    // x = x + K * y
    for (int i = 0; i < STATE_DIM; i++)
    {
        float temp = 0;
        for (int j = 0; j < MEASURE_DIM; j++)
        {
            temp += K[i][j] * y[j];
        }
        kf->x_state[i] += temp;
    }

    // P = (I - K * H) * P
    // Costruisci I - K * H e poi moltiplica con P
    float I[STATE_DIM][STATE_DIM];
    for (int i = 0; i < STATE_DIM; i++)
    {
        for (int j = 0; j < STATE_DIM; j++)
        {
            I[i][j] = (i == j) ? 1 : 0;
        }
    }
    // matrix_multiply(K, H, temp1);
    // matrix_subtract(I, temp1, temp1);
    // matrix_multiply(temp1, P, P);
    matrix_multiply(STATE_DIM, MEASURE_DIM, STATE_DIM, K, kf->H, temp1);
    matrix_subtract(STATE_DIM, STATE_DIM, I, temp1, temp1);
    matrix_multiply(STATE_DIM, STATE_DIM, STATE_DIM, temp1, kf->P, kf->P);

    // print x,y,z
    // printf("x = %f, y = %f, z = %f\n", kf->x_state[0], kf->x_state[1], kf->x_state[2]);
}