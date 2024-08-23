#ifndef LINEAR_KALMAN_FILTER_RW_H
#define LINEAR_KALMAN_FILTER_RW_H

#include "arm_math.h" // Libreria CMSIS-DSP per operazioni con matrici
#include "MagneticDeck.h"

// Dichiarazione delle matrici del Filtro di Kalman
extern float32_t F_data[STATE_DIM * STATE_DIM];
extern float32_t H_data[MEASURE_DIM * STATE_DIM];
extern float32_t Q_data[STATE_DIM * STATE_DIM];
extern float32_t R_data[MEASURE_DIM * MEASURE_DIM];
extern float32_t P_data[STATE_DIM * STATE_DIM];
extern float32_t x_data[STATE_DIM];
extern float32_t z_data[MEASURE_DIM];

// Dichiarazione degli oggetti matrice CMSIS-DSP
extern arm_matrix_instance_f32 F, H, Q, R, P, x, z, K;

// Dichiarazione delle funzioni

void kalman_init(KalmanFilter *kf, float *initial_guess);
void kalman_predict(KalmanFilter *kf);
void kalman_update(KalmanFilter *kf);

#endif // LINEAR_KALMAN_FILTER_RW_H