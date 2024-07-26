#include "math.h"
#include "mm_magnetic_distance.h"
#include "log.h"
#include "debug.h"
#include "param.h"
#include "stabilizer_types.h"
#include "estimator.h"
#include "estimator_kalman.h"
#include <time.h>
#include <stdlib.h>
#include <time.h>
#include <MagneticDeck.h>

#include "usec_time.h"

#include "nelder_mead_4A.h"
#include "nelder_mead_3A.h"

float Optimization_Model_STD = 0.05f;

// #define start_timer() *((volatile uint32_t *)0xE0001000) = 0x40000001 // Enable CYCCNT register
// #define stop_timer() *((volatile uint32_t *)0xE0001000) = 0x40000000  // Disable CYCCNT register
// #define get_timer() *((volatile uint32_t *)0xE0001004)

// logVarId_t gainID =  logGetVarId("Potentiometer_G_P", "GainValue");
// #define G_INA logGetFloat(GainID)

// uint32_t it1, it2; // start and stop flag
volatile float MeasuredVoltages_calibrated[4] = {0.0f, 0.0f, 0.0f, 0.0f};
float PredictedVoltages[4] = {0.0f, 0.0f, 0.0f, 0.0f};

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
volatile float InputPoint[3] = {0.0f, 0.0f, 0.0f};
volatile float noise = 0.0f;

float default_CG_a1 = 1.3f; //  NERO
float CG_a1 = 1.0f;         //  NERO

float default_CG_a2 = 0.635f; // GIALLO
float CG_a2 = 1.0f;           // GIALLO

float default_CG_a3 = 3.5f; // GRIGIO
float CG_a3 = 1.0f;         // GRIGIO

float default_CG_a4 = 7.5f; // ROSSO
float CG_a4 = 1.0f;         // ROSSO

float calibrationsGains[4];

float CalibrationRapport_anchor1 = 1.0f;
float CalibrationRapport_anchor2 = 1.0f;
float CalibrationRapport_anchor3 = 1.0f;
float CalibrationRapport_anchor4 = 1.0f;

// ----------------- Calibration -----------------
volatile uint32_t currentCalibrationTick = 0;
float calibrationMean[4] = {};

// adaptive std on Measured volts
int number_of_samples_for_DynamicStd = 0;
float Nero_std_data[window_size];
float Giallo_std_data[window_size];
float Grigio_std_data[window_size];
float Rosso_std_data[window_size];

volatile float Nero_std = 0;
volatile float Giallo_std = 0;
volatile float Grigio_std = 0;
volatile float Rosso_std = 0;

static float B_field_vector_1[3];
static float B_field_vector_2[3];
static float B_field_vector_3[3];
static float B_field_vector_4[3];

static bool isInitialized = false;
static bool isFirstMeasurement = true;

static float estimated_position[3];

static int counter = 0;

// Set the range where to look for the minumum
static float range[3] = {-0.2f, +0.2f, 0.2f};
static float solution[3];
static nm_params_t_3A paramsNM_3A;
static nm_params_t_4A paramsNM_4A;

void kalmanCoreUpdateWithVolt(kalmanCoreData_t *this, voltMeasurement_t *voltAnchor)
{

    // float all_start = usecTimestamp();

    if (!isInitialized)
    {
        // initialize the random seed
        srand(time(NULL));
        isInitialized = true;

        // initialize the optimization model parameters

        nm_params_init_default_3A(&paramsNM_3A, 3);
        paramsNM_3A.debug_log = 0;
        paramsNM_3A.max_iterations = 15;
        paramsNM_3A.tol_fx = 1e-6f;
        paramsNM_3A.tol_x = 1e-5f;
        paramsNM_3A.restarts = 0;

        nm_params_init_default_4A(&paramsNM_4A, 3);
        paramsNM_4A.debug_log = 0;
        paramsNM_4A.max_iterations = 15;
        paramsNM_4A.tol_fx = 1e-6f;
        paramsNM_4A.tol_x = 1e-5f;
        paramsNM_4A.restarts = 0;
    }

    if (currentCalibrationTick < CALIBRATION_TIC_VALUE)
    {
        // ------------------------------ COLLECTING CALIBRATION MEASUREMENTS ------------------------------
        calibrationMean[0] += voltAnchor->measuredVolt[0];
        calibrationMean[1] += voltAnchor->measuredVolt[1];
        calibrationMean[2] += voltAnchor->measuredVolt[2];
        calibrationMean[3] += voltAnchor->measuredVolt[3];

        currentCalibrationTick = currentCalibrationTick + 1;
    }
    else
    {
        if (currentCalibrationTick == CALIBRATION_TIC_VALUE)
        {
            // ------------------------------ CALIBRATION ------------------------------
            // For  calibration i fixed the tag position and orientation to 0,0,0 and 0,0,1
            currentCalibrationTick = currentCalibrationTick + 1;

            float meanData_a1 = calibrationMean[0] / CALIBRATION_TIC_VALUE;
            float meanData_a2 = calibrationMean[1] / CALIBRATION_TIC_VALUE;
            float meanData_a3 = calibrationMean[2] / CALIBRATION_TIC_VALUE;
            float meanData_a4 = calibrationMean[3] / CALIBRATION_TIC_VALUE;

            // DEBUG_PRINT("calibration data acquired!!\n");
            DEBUG_PRINT("calibration_Mean[0][0] = %f\n", (double)meanData_a1);
            DEBUG_PRINT("calibration_Mean[1][1] = %f\n", (double)meanData_a2);
            DEBUG_PRINT("calibration_Mean[2][2] = %f\n", (double)meanData_a3);
            DEBUG_PRINT("calibration_Mean[3][3] = %f\n", (double)meanData_a4);

            // point_t cfPosP;
            // estimatorKalmanGetEstimatedPos(&cfPosP);
            // float tag_pos_predicted_calibrated[3] = {cfPosP.x, cfPosP.y, cfPosP.z};
            float tag_pos_predicted_calibrated[3] = {0.0f, 0.0f, 0.0f};

            // float RotationMatrix[3][3];
            // estimatorKalmanGetEstimatedRot((float *)RotationMatrix);
            // float tag_or_versor_calibrated[3] = {RotationMatrix[0][2], RotationMatrix[1][2], RotationMatrix[2][2]};
            float tag_or_versor_calibrated[3] = {0.0f, 0.0f, 1.0f};

            float anchor_1_pose[3] = {voltAnchor->x[0], voltAnchor->y[0], voltAnchor->z[0]};
            float anchor_2_pose[3] = {voltAnchor->x[1], voltAnchor->y[1], voltAnchor->z[1]};
            float anchor_3_pose[3] = {voltAnchor->x[2], voltAnchor->y[2], voltAnchor->z[2]};
            float anchor_4_pose[3] = {voltAnchor->x[3], voltAnchor->y[3], voltAnchor->z[3]};
            get_B_field_for_a_Anchor(anchor_1_pose, tag_pos_predicted_calibrated, tag_or_versor_calibrated, B_field_vector_1);
            get_B_field_for_a_Anchor(anchor_2_pose, tag_pos_predicted_calibrated, tag_or_versor_calibrated, B_field_vector_2);
            get_B_field_for_a_Anchor(anchor_3_pose, tag_pos_predicted_calibrated, tag_or_versor_calibrated, B_field_vector_3);
            get_B_field_for_a_Anchor(anchor_4_pose, tag_pos_predicted_calibrated, tag_or_versor_calibrated, B_field_vector_4);

            // computing the V_rx for each of the 4 anchors
            float V_rx_1 = V_from_B(B_field_vector_1, tag_or_versor_calibrated, voltAnchor->resonanceFrequency[0], voltAnchor->GainValue);
            float V_rx_2 = V_from_B(B_field_vector_2, tag_or_versor_calibrated, voltAnchor->resonanceFrequency[1], voltAnchor->GainValue);
            float V_rx_3 = V_from_B(B_field_vector_3, tag_or_versor_calibrated, voltAnchor->resonanceFrequency[2], voltAnchor->GainValue);
            float V_rx_4 = V_from_B(B_field_vector_4, tag_or_versor_calibrated, voltAnchor->resonanceFrequency[3], voltAnchor->GainValue);

            CG_a1 = meanData_a1 / V_rx_1;
            CG_a2 = meanData_a2 / V_rx_2;
            CG_a3 = meanData_a3 / V_rx_3;
            CG_a4 = meanData_a4 / V_rx_4;
            calibrationsGains[0] = CG_a1;
            calibrationsGains[1] = CG_a2;
            calibrationsGains[2] = CG_a3;
            calibrationsGains[3] = CG_a4;

            DEBUG_PRINT("CG_a1 = %f\n", (double)CG_a1);
            DEBUG_PRINT("CG_a2 = %f\n", (double)CG_a2);
            DEBUG_PRINT("CG_a3 = %f\n", (double)CG_a3);
            DEBUG_PRINT("CG_a4 = %f\n", (double)CG_a4);
        }
        else
        {

            if (currentCalibrationTick == CALIBRATION_TIC_VALUE + 1)
            {
                DEBUG_PRINT("currentCalibrationTick = %f\n", currentCalibrationTick);

                DEBUG_PRINT("Resetting the Kalman filte after calibrationr\n");
                paramSetInt(paramGetVarId("kalman", "resetEstimation"), 1);
                currentCalibrationTick = currentCalibrationTick + 1;
                paramSetInt(paramGetVarId("kalman", "resetEstimation"), 0);

                // paramSetInt(paramGetVarId("kalman", "resetEstimation"), 0);
            }

            // ------------------------------------ RUNNING CASE ------------------------------------

            float RotationMatrix[3][3];
            // it would be the product between the rotation matrix and the initial [0,0,1] versor
            estimatorKalmanGetEstimatedRot((float *)RotationMatrix);
            float tag_or_versor[3] = {RotationMatrix[0][2], RotationMatrix[1][2], RotationMatrix[2][2]};

            // for (int anchorIdx = 0; anchorIdx < NUM_ANCHORS; anchorIdx++)
            // {
            // }

            /// ------------------------- optimization  for computing position -------------------------

            // initialize the params for the optimization algorithm
            // float my_params_start_cost_all = usecTimestamp();
            float x_start[3] = {};
            if (isFirstMeasurement)
            {
                x_start[0] = this->S[KC_STATE_X];
                x_start[1] = this->S[KC_STATE_Y];
                x_start[2] = this->S[KC_STATE_Z];
                isFirstMeasurement = false;
            }
            else
            {
                // TODO: forse qui si puo fare di meglio che usare la stima precedente.
                // Si potrebbe ad esempio controllare se la stima precedente è un outlier e nel caso scartarla e usare lo stato de kalman come stima iniziale,
                // o una predizione basata su qualche modello lineare di n stime precedenti
                // QUESTI CONTROLLI VANNO FATTI PRIMA DI FARE UN UPDATE DEL KALMAN
                x_start[0] = estimated_position[0];
                x_start[1] = estimated_position[1];
                x_start[2] = estimated_position[2];
            }

            if (voltAnchor->Id_in_saturation != 10)
            {
                myParams_t_3A my_params;
                for (int anchorIdx = 0; anchorIdx < NUM_ANCHORS; anchorIdx++)
                {
                    my_params.Anchors_3A[anchorIdx][0] = voltAnchor->x[anchorIdx];
                    my_params.Anchors_3A[anchorIdx][1] = voltAnchor->y[anchorIdx];
                    my_params.Anchors_3A[anchorIdx][2] = voltAnchor->z[anchorIdx];

                    my_params.versore_orientamento_cf_3A[0] = tag_or_versor[0];
                    my_params.versore_orientamento_cf_3A[1] = tag_or_versor[1];
                    my_params.versore_orientamento_cf_3A[2] = tag_or_versor[2];

                    my_params.Gain_3A = voltAnchor->GainValue;

                    my_params.frequencies_3A[anchorIdx] = voltAnchor->resonanceFrequency[anchorIdx];
                    // using the measured volts to fill the structure
                    my_params.MeasuredVoltages_calibrated_3A[anchorIdx] = voltAnchor->measuredVolt[anchorIdx] / calibrationsGains[anchorIdx];
                    MeasuredVoltages_calibrated[anchorIdx] = my_params.MeasuredVoltages_calibrated_3A[anchorIdx];
                }
                // float my_params_ms = (float)(usecTimestamp() - my_params_start_cost_all) / 1000.0f;

                // set the starting point for the optimization algorithm

                // float optimize_start_cost_all = usecTimestamp();
                nm_result_t_3A result = nm_multivar_optimize_3A(3, x_start, range, &myCostFunction_3A, &my_params, &paramsNM_3A, solution);
                // float optimize_ms = (float)(usecTimestamp() - optimize_start_cost_all) / 1000.0f;
            }
            else
            {
                myParams_t_4A my_params;
                for (int anchorIdx = 0; anchorIdx < NUM_ANCHORS; anchorIdx++)
                {
                    my_params.Anchors_4A[anchorIdx][0] = voltAnchor->x[anchorIdx];
                    my_params.Anchors_4A[anchorIdx][1] = voltAnchor->y[anchorIdx];
                    my_params.Anchors_4A[anchorIdx][2] = voltAnchor->z[anchorIdx];

                    my_params.versore_orientamento_cf_4A[0] = tag_or_versor[0];
                    my_params.versore_orientamento_cf_4A[1] = tag_or_versor[1];
                    my_params.versore_orientamento_cf_4A[2] = tag_or_versor[2];

                    my_params.Gain_4A = voltAnchor->GainValue;

                    my_params.frequencies_4A[anchorIdx] = voltAnchor->resonanceFrequency[anchorIdx];
                    // using the measured volts to fill the structure
                    my_params.MeasuredVoltages_calibrated_4A[anchorIdx] = voltAnchor->measuredVolt[anchorIdx] / calibrationsGains[anchorIdx];
                    MeasuredVoltages_calibrated[anchorIdx] = my_params.MeasuredVoltages_calibrated_4A[anchorIdx];
                }
                // float my_params_ms = (float)(usecTimestamp() - my_params_start_cost_all) / 1000.0f;

                // set the starting point for the optimization algorithm

                // float optimize_start_cost_all = usecTimestamp();
                nm_result_t_4A result = nm_multivar_optimize_4A(3, x_start, range, &myCostFunction_4A, &my_params, &paramsNM_4A, solution);
                // float optimize_ms = (float)(usecTimestamp() - optimize_start_cost_all) / 1000.0f;
            }

            /// ------------------------- optimization  for computing position -------------------------

            float h_x[KC_STATE_DIM] = {0};
            arm_matrix_instance_f32 H_x = {1, KC_STATE_DIM, h_x};
            h_x[KC_STATE_X] = 1;

            float h_y[KC_STATE_DIM] = {0};
            arm_matrix_instance_f32 H_y = {1, KC_STATE_DIM, h_y};
            h_y[KC_STATE_Y] = 1;

            float h_z[KC_STATE_DIM] = {0};
            arm_matrix_instance_f32 H_z = {1, KC_STATE_DIM, h_z};
            h_z[KC_STATE_Z] = 1;

            // float EKF_UPDATE_start_cost_all = usecTimestamp();
            kalmanCoreScalarUpdate(this, &H_x, solution[0] - this->S[KC_STATE_X], Optimization_Model_STD);
            kalmanCoreScalarUpdate(this, &H_y, solution[1] - this->S[KC_STATE_Y], Optimization_Model_STD);
            kalmanCoreScalarUpdate(this, &H_z, solution[2] - this->S[KC_STATE_Z], Optimization_Model_STD);
            // float EKF_UPDATE_ms = (float)(usecTimestamp() - EKF_UPDATE_start_cost_all) / 1000.0f;

            // float ALL_ms = (float)(usecTimestamp() - all_start) / 1000.0f;
            // DEBUG_PRINT("ALL_ms = %f\n", (double)ALL_ms);

            // float a = 0;
            // DEBUG_PRINT("Estimated position x: %f,\n", (double)estimated_position[0]);
        }
    }
}

LOG_GROUP_START(Optimization_Model)
LOG_ADD(LOG_FLOAT, T_x, &solution[0])
LOG_ADD(LOG_FLOAT, T_y, &solution[1])
LOG_ADD(LOG_FLOAT, T_z, &solution[2])

// LOG_ADD(LOG_FLOAT, Inpt_x, &InputPoint[0])
// LOG_ADD(LOG_FLOAT, Inpt_y, &InputPoint[1])
// LOG_ADD(LOG_FLOAT, Inpt_z, &InputPoint[2])

LOG_GROUP_STOP(Optimization_Model)

// PARAM_GROUP_START(Opt_Model_Param)
// PARAM_ADD(PARAM_FLOAT, Opt_STD, &Optimization_Model_STD)
// PARAM_GROUP_STOP(Opt_Model_Param)

LOG_GROUP_START(Dipole_Model)

// // LOG_ADD(LOG_UINT32, CPUCycle, &it2)

// LOG_ADD(LOG_FLOAT, CG_A1, &CG_a1)
// LOG_ADD(LOG_FLOAT, CG_A2, &CG_a2)
// LOG_ADD(LOG_FLOAT, CG_A3, &CG_a3)
// LOG_ADD(LOG_FLOAT, CG_A4, &CG_a4)

// LOG_ADD(LOG_FLOAT, ERR1, &E_1)
// LOG_ADD(LOG_FLOAT, ERR2, &E_2)
// LOG_ADD(LOG_FLOAT, ERR3, &E_3)
// LOG_ADD(LOG_FLOAT, ERR4, &E_4)

LOG_ADD(LOG_FLOAT, M_V1, &MeasuredVoltages_calibrated[0])
LOG_ADD(LOG_FLOAT, M_V2, &MeasuredVoltages_calibrated[1])
LOG_ADD(LOG_FLOAT, M_V3, &MeasuredVoltages_calibrated[2])
LOG_ADD(LOG_FLOAT, M_V4, &MeasuredVoltages_calibrated[3])

// LOG_ADD(LOG_FLOAT, P_V_0, &PredictedVoltages[0])
// LOG_ADD(LOG_FLOAT, P_V_1, &PredictedVoltages[1])
// LOG_ADD(LOG_FLOAT, P_V_2, &PredictedVoltages[2])
// LOG_ADD(LOG_FLOAT, P_V_3, &PredictedVoltages[3])

LOG_GROUP_STOP(Dipole_Model)

// PARAM_GROUP_START(Dipole_Params)
// PARAM_ADD(PARAM_UINT16, calibTic, &currentCalibrationTick)
// PARAM_GROUP_STOP(Dipole_Params)