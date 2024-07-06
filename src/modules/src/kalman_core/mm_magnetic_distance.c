#include "math.h"
#include "mm_magnetic_distance.h"
#include "log.h"
#include "debug.h"
#include "param.h"
#include "stabilizer_types.h"
#include "estimator.h"
#include "estimator_kalman.h"
#include "nelder_mead.h"

#define G_INA 100.0f
#define RAY 0.019f
#define N_WOUNDS 5.0f
#define COIL_SURFACE (RAY * RAY * PI)
#define CURRENT 0.5f
#define MU_0 1.25663706212e-06f
#define PI 3.14159265359f
#define Num_Anchors 4

// #define start_timer() *((volatile uint32_t *)0xE0001000) = 0x40000001 // Enable CYCCNT register
// #define stop_timer() *((volatile uint32_t *)0xE0001000) = 0x40000000  // Disable CYCCNT register
// #define get_timer() *((volatile uint32_t *)0xE0001004)

// logVarId_t gainID =  logGetVarId("Potentiometer_G_P", "GainValue");
// #define G_INA logGetFloat(GainID)

// uint32_t it1, it2; // start and stop flag
volatile float MeasuredVoltages_calibrated[4] = {0.0f, 0.0f, 0.0f, 0.0f};
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

float T_pre_x_fk = 0.0f;
float T_pre_y_fk = 0.0f;
float T_pre_z_fk = 0.0f;

float T_orientation[3] = {0.0f, 0.0f, 0.0f};

float default_CG_a1 = 1.3f; //  NERO
float CG_a1 = 1.0f;         //  NERO

float default_CG_a2 = 0.635f; // GIALLO
float CG_a2 = 1.0f;           // GIALLO

float default_CG_a3 = 3.5f; // GRIGIO
float CG_a3 = 1.0f;         // GRIGIO

float default_CG_a4 = 7.5f; // ROSSO
float CG_a4 = 1.0f;         // ROSSO

float CalibrationRapport_anchor1 = 1.0f;
float CalibrationRapport_anchor2 = 1.0f;
float CalibrationRapport_anchor3 = 1.0f;
float CalibrationRapport_anchor4 = 1.0f;

volatile uint32_t pos_idx = 0;

// ----------------- Calibration -----------------
#define CALIBRATION_TIC_VALUE 100.0f // Replace with the desired constant value
volatile uint32_t currentCalibrationTick = 0;
static float calibrationMean[4][4] = {};

// adaptive std on Measured volts
#define window_size 25
int number_of_samples_for_DynamicStd = 0;
float Nero_std_data[window_size];
float Giallo_std_data[window_size];
float Grigio_std_data[window_size];
float Rosso_std_data[window_size];

volatile float Nero_std = 0;
volatile float Giallo_std = 0;
volatile float Grigio_std = 0;
volatile float Rosso_std = 0;

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

optimset opt = {
    .precision = 5.0f,
    .format = 0.0f,
    .verbose = 0.0f,
    .tol_x = 1e-3f,
    .tol_y = 1e-3f,
    .max_iter = 500.0f,
    .max_eval = 500.0f,
    .adaptive = 0.0f,
    .scale = 1.0e-3f};

const int problem_dimension = 3;

struct Model
{
    real V[Num_Anchors];
};

void init_model(float *tag_pos_predicted, float *tag_or_versor, voltMeasurement_t *voltAnchor, model *mdl)
{
    // model *mdl = malloc(Num_Anchors * sizeof(model));
    ///============================= initialization of the stuff for the magnetic simulation ==============================

    // calculate the B field and the V for each coil in each positions
    // this is the ground truth
    real B_field[3];
    real V[Num_Anchors];

    for (int anchorIdx = 0; anchorIdx < Num_Anchors; anchorIdx++)
    {
        real Anchor[3] = {voltAnchor->x[anchorIdx], voltAnchor->y[anchorIdx], voltAnchor->z[anchorIdx]};
        get_B_field_for_a_Anchor(Anchor, tag_pos_predicted, tag_or_versor, B_field);
        V[anchorIdx] = V_from_B(B_field, tag_or_versor, voltAnchor->resonanceFrequency[anchorIdx]);
        mdl->V[anchorIdx] = V[anchorIdx];
    }
}

int dimensions()
{
    return problem_dimension; // Tre componenti per la posizione (x, y, z)
}

void cost(const model *mdl, point *pnt)
{
    real costo = 0.0f;
    for (int i = 0; i < Num_Anchors; i++)
    {
        costo += powf(mdl->V[i] - MeasuredVoltages_calibrated[i], 2);
    }

    pnt->y = costo;
}

float estimated_position[3];
float B_field_vector_1[3];
float B_field_vector_2[3];
float B_field_vector_3[3];
float B_field_vector_4[3];

point inp;
point out;
model mdl;
simplex smpl;

void kalmanCoreUpdateWithVolt(kalmanCoreData_t *this, voltMeasurement_t *voltAnchor)
{

    if (currentCalibrationTick < CALIBRATION_TIC_VALUE)
    {
        // ------------------------------ COLLECTING CALIBRATION MEASUREMENTS ------------------------------
        calibrationMean[0][0] += voltAnchor->measuredVolt[0];
        calibrationMean[1][1] += voltAnchor->measuredVolt[1];
        calibrationMean[2][2] += voltAnchor->measuredVolt[2];
        calibrationMean[3][3] += voltAnchor->measuredVolt[3];

        currentCalibrationTick = currentCalibrationTick + 1;
    }
    else
    {
        if (currentCalibrationTick == CALIBRATION_TIC_VALUE)
        {
            // ------------------------------ CALIBRATION ------------------------------
            // For  calibration i fixed the tag position and orientation to 0,0,0 and 0,0,1
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

            point_t cfPosP;
            estimatorKalmanGetEstimatedPos(&cfPosP);
            // float tag_pos_predicted_calibrated[3] = {cfPosP.x, cfPosP.y, cfPosP.z};
            float tag_pos_predicted_calibrated[3] = {0.0, 0.0, 0.0};
            T_pre_x_fk = tag_pos_predicted_calibrated[0];
            T_pre_y_fk = tag_pos_predicted_calibrated[1];
            T_pre_z_fk = tag_pos_predicted_calibrated[2];

            float RotationMatrix[3][3];
            estimatorKalmanGetEstimatedRot((float *)RotationMatrix);
            // float tag_or_versor_calibrated[3] = {RotationMatrix[0][2], RotationMatrix[1][2], RotationMatrix[2][2]};
            float tag_or_versor_calibrated[3] = {0.0, 0.0, 1.0};

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
            // DEBUG_PRINT("sRunning case!!\n");
            // ------------------------------------ RUNNING CASE ------------------------------------
            // start_timer();     // start the timer.
            // it1 = get_timer(); // store current cycle-count in a local

            // computing the B field for each of the 4 anchors
            point_t cfPosP;
            estimatorKalmanGetEstimatedPos(&cfPosP);
            float cfPos[3] = {cfPosP.x, cfPosP.y, cfPosP.z};
            T_pre_x_fk = cfPos[0];
            T_pre_y_fk = cfPos[1];
            T_pre_z_fk = cfPos[2];

            float tag_pos_predicted[3] = {cfPos[0], cfPos[1], cfPos[2]};
            // float tag_pos_predicted[3] = {0.0, 0.0, 0.0};

            float RotationMatrix[3][3];
            // it would be the product between the rotation matrix and the initial [0,0,1] versor
            estimatorKalmanGetEstimatedRot((float *)RotationMatrix);
            float tag_or_versor[3] = {RotationMatrix[0][2], RotationMatrix[1][2], RotationMatrix[2][2]};
            // float tag_or_versor[3] = {0.0, 0.0, 1.0};

            // T_orientation[0] = tag_or_versor[0];
            // T_orientation[1] = tag_or_versor[1];
            // T_orientation[2] = tag_or_versor[2];

            // float anchor_1_pose[3] = {voltAnchor->x[0], voltAnchor->y[0], voltAnchor->z[0]};
            // float anchor_2_pose[3] = {voltAnchor->x[1], voltAnchor->y[1], voltAnchor->z[1]};
            // float anchor_3_pose[3] = {voltAnchor->x[2], voltAnchor->y[2], voltAnchor->z[2]};
            // float anchor_4_pose[3] = {voltAnchor->x[3], voltAnchor->y[3], voltAnchor->z[3]};

            // MeasuredVoltages_raw[0] = voltAnchor->measuredVolt[0];
            // MeasuredVoltages_raw[1] = voltAnchor->measuredVolt[1];
            // MeasuredVoltages_raw[2] = voltAnchor->measuredVolt[2];
            // MeasuredVoltages_raw[3] = voltAnchor->measuredVolt[3];

            MeasuredVoltages_calibrated[0] = (voltAnchor->measuredVolt[0] / CG_a1) / 0.995941558f;
            MeasuredVoltages_calibrated[1] = voltAnchor->measuredVolt[1] / CG_a2;
            MeasuredVoltages_calibrated[2] = voltAnchor->measuredVolt[2] / CG_a3;
            MeasuredVoltages_calibrated[3] = voltAnchor->measuredVolt[3] / CG_a4;

            /// ------------------------- optimization  for computing position -------------------------

            const int n = problem_dimension;

            // initialize the initial optimization point
            inp.x[0] = cfPos[0];
            inp.x[1] = cfPos[1];
            inp.x[2] = cfPos[2];

            // initialize model and simplex

            init_model(tag_pos_predicted, tag_or_versor, voltAnchor, &mdl);

            init_simplex(n, opt.scale, &inp, &smpl);

            cost(&mdl, &inp);

            // optimize model mdl, using settings opt, starting
            // from simplex smpl, and save result to point out
            nelder_mead(&mdl, &opt, &smpl, &out);

            estimated_position[0] = out.x[0];
            estimated_position[1] = out.x[1];
            estimated_position[2] = out.x[2];
            // free_simplex(&smpl);
            // free_point(&inp);
            // free_point(&out);
            // free(&mdl);

            /// ------------------------- optimization  for computing position -------------------------

            // float error_x = ;
            // float error_y = ;
            // float error_z = ;

            float h_x[KC_STATE_DIM] = {0};
            arm_matrix_instance_f32 H_x = {1, KC_STATE_DIM, h_x};
            h_x[KC_STATE_X] = 1;

            float h_y[KC_STATE_DIM] = {0};
            arm_matrix_instance_f32 H_y = {1, KC_STATE_DIM, h_y};
            h_y[KC_STATE_Y] = 1;

            float h_z[KC_STATE_DIM] = {0};
            arm_matrix_instance_f32 H_z = {1, KC_STATE_DIM, h_z};
            h_z[KC_STATE_Z] = 1;

            kalmanCoreScalarUpdate(this, &H_x, estimated_position[0] - cfPos[0], 0.05f);
            kalmanCoreScalarUpdate(this, &H_y, estimated_position[1] - cfPos[1], 0.05f);
            kalmanCoreScalarUpdate(this, &H_z, estimated_position[2] - cfPos[2], 0.05f);

            // DEBUG_PRINT("Estimated position x: %f,\n", (double)estimated_position[0]);
        }
    }
}

LOG_GROUP_START(Optimization_Model)
LOG_ADD(LOG_FLOAT, T_x, &estimated_position[0])
LOG_ADD(LOG_FLOAT, T_y, &estimated_position[1])
LOG_ADD(LOG_FLOAT, T_z, &estimated_position[2])
LOG_GROUP_STOP(Optimization_Model)

// LOG_GROUP_START(Dipole_Model)

// // LOG_ADD(LOG_UINT32, CPUCycle, &it2)

// LOG_ADD(LOG_FLOAT, D_x_0, &Derivative_x[0])
// LOG_ADD(LOG_FLOAT, D_x_1, &Derivative_x[1])
// LOG_ADD(LOG_FLOAT, D_x_2, &Derivative_x[2])
// LOG_ADD(LOG_FLOAT, D_x_3, &Derivative_x[3])
// LOG_ADD(LOG_FLOAT, D_y_0, &Derivative_y[0])
// LOG_ADD(LOG_FLOAT, D_y_1, &Derivative_y[1])
// LOG_ADD(LOG_FLOAT, D_y_2, &Derivative_y[2])
// LOG_ADD(LOG_FLOAT, D_y_3, &Derivative_y[3])
// LOG_ADD(LOG_FLOAT, D_z_0, &Derivative_z[0])
// LOG_ADD(LOG_FLOAT, D_z_1, &Derivative_z[1])
// LOG_ADD(LOG_FLOAT, D_z_2, &Derivative_z[2])
// LOG_ADD(LOG_FLOAT, D_z_3, &Derivative_z[3])

// LOG_ADD(LOG_FLOAT, CG_A1, &CG_a1)
// LOG_ADD(LOG_FLOAT, CG_A2, &CG_a2)
// LOG_ADD(LOG_FLOAT, CG_A3, &CG_a3)
// LOG_ADD(LOG_FLOAT, CG_A4, &CG_a4)

// LOG_ADD(LOG_FLOAT, C_R_1, &CalibrationRapport_anchor1)
// LOG_ADD(LOG_FLOAT, C_R_2, &CalibrationRapport_anchor2)
// LOG_ADD(LOG_FLOAT, C_R_3, &CalibrationRapport_anchor3)
// LOG_ADD(LOG_FLOAT, C_R_4, &CalibrationRapport_anchor4)

// // LOG_ADD(LOG_FLOAT, T_pred_x, &T_pre_x)
// // LOG_ADD(LOG_FLOAT, T_pred_y, &T_pre_y)
// // LOG_ADD(LOG_FLOAT, T_pred_z, &T_pre_z)

// LOG_ADD(LOG_FLOAT, T_pred_x_fk, &T_pre_x_fk)
// LOG_ADD(LOG_FLOAT, T_pred_y_fk, &T_pre_y_fk)
// LOG_ADD(LOG_FLOAT, T_pred_z_fk, &T_pre_z_fk)

// LOG_ADD(LOG_FLOAT, T_or_0, &T_orientation[0])
// LOG_ADD(LOG_FLOAT, T_or_1, &T_orientation[1])
// LOG_ADD(LOG_FLOAT, T_or_2, &T_orientation[2])

// LOG_ADD(LOG_FLOAT, Er_1, &E_1)
// LOG_ADD(LOG_FLOAT, Er_2, &E_2)
// LOG_ADD(LOG_FLOAT, Er_3, &E_3)
// LOG_ADD(LOG_FLOAT, Er_4, &E_4)

// LOG_ADD(LOG_FLOAT, MC_V_0, &MeasuredVoltages_calibrated[0])
// LOG_ADD(LOG_FLOAT, MC_V_1, &MeasuredVoltages_calibrated[1])
// LOG_ADD(LOG_FLOAT, MC_V_2, &MeasuredVoltages_calibrated[2])
// LOG_ADD(LOG_FLOAT, MC_V_3, &MeasuredVoltages_calibrated[3])

// LOG_ADD(LOG_FLOAT, MR_V_0, &MeasuredVoltages_raw[0])
// LOG_ADD(LOG_FLOAT, MR_V_1, &MeasuredVoltages_raw[1])
// LOG_ADD(LOG_FLOAT, MR_V_2, &MeasuredVoltages_raw[2])
// LOG_ADD(LOG_FLOAT, MR_V_3, &MeasuredVoltages_raw[3])

// LOG_ADD(LOG_FLOAT, P_V_0, &PredictedVoltages[0])
// LOG_ADD(LOG_FLOAT, P_V_1, &PredictedVoltages[1])
// LOG_ADD(LOG_FLOAT, P_V_2, &PredictedVoltages[2])
// LOG_ADD(LOG_FLOAT, P_V_3, &PredictedVoltages[3])

// LOG_ADD(LOG_FLOAT, B_0_0, &B[0][0])
// LOG_ADD(LOG_FLOAT, B_0_1, &B[0][1])
// LOG_ADD(LOG_FLOAT, B_0_2, &B[0][2])

// LOG_ADD(LOG_FLOAT, B_1_0, &B[1][0])
// LOG_ADD(LOG_FLOAT, B_1_1, &B[1][1])
// LOG_ADD(LOG_FLOAT, B_1_2, &B[1][2])

// LOG_ADD(LOG_FLOAT, B_2_0, &B[2][0])
// LOG_ADD(LOG_FLOAT, B_2_1, &B[2][1])
// LOG_ADD(LOG_FLOAT, B_2_2, &B[2][2])

// LOG_ADD(LOG_FLOAT, B_3_0, &B[3][0])
// LOG_ADD(LOG_FLOAT, B_3_1, &B[3][1])
// LOG_ADD(LOG_FLOAT, B_3_2, &B[3][2])

// LOG_ADD(LOG_FLOAT, A0_N_std, &Nero_std)
// LOG_ADD(LOG_FLOAT, A1_Gia_std, &Giallo_std)
// LOG_ADD(LOG_FLOAT, A2_Gr_std, &Grigio_std)
// LOG_ADD(LOG_FLOAT, A3_R_std, &Rosso_std)

// LOG_ADD(LOG_FLOAT, A3_R_std, &Rosso_std)

// LOG_GROUP_STOP(Dipole_Model)

// PARAM_GROUP_START(Dipole_Params)
// PARAM_ADD(PARAM_UINT16, calibTic, &currentCalibrationTick)
// PARAM_GROUP_STOP(Dipole_Params)