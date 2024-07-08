typedef struct
{
    float P[3][3];
    float x[3];
    float R;
    float Q;
} ekf_t;

void kalman_custom_update(ekf_t *ekf, float H[3], float delta);
void kalman_custom_predict(ekf_t *ekf);