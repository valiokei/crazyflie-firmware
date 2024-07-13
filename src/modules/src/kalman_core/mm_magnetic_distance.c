#include "math.h"
#include "mm_magnetic_distance.h"
#include "log.h"
#include "debug.h"
#include "param.h"
#include "stabilizer_types.h"
#include "estimator.h"
#include "estimator_kalman.h"
#include <stdlib.h>
#include "MagneticDeck.h"

// ------------------------------ Variables for computing execution time ------------------------------
// #define start_timer() *((volatile uint32_t *)0xE0001000) = 0x40000001 // Enable CYCCNT register
// #define stop_timer() *((volatile uint32_t *)0xE0001000) = 0x40000000  // Disable CYCCNT register
// #define get_timer() *((volatile uint32_t *)0xE0001004)

// uint32_t it1, it2; // start and stop flag

// ------------------------------ Calibration ------------------------------
float CG_a1 = 1.0f; //  NERO
float CG_a2 = 1.0f; // GIALLO
float CG_a3 = 1.0f; // GRIGIO
float CG_a4 = 1.0f; // ROSSO
float currentCalibrationTick = 0.0;
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

// ------------------------------ Dipole Model  Variables Initialization------------------------------
float B_field_vector_1[3];
float B_field_vector_2[3];
float B_field_vector_3[3];
float B_field_vector_4[3];

float V_rx_derivative_x[4];
float V_rx_derivative_y[4];
float V_rx_derivative_z[4];

// ------------------------------ DEBUG ------------------------------
static float PredictedVoltages[4];
static float MeasuredVoltages[4];
static int printFlagSTD_Type = 0;

void kalmanCoreUpdateWithVolt(kalmanCoreData_t *this, voltMeasurement_t *voltAnchor)
{

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
            float tag_pos_predicted_calibrated[3] = {0.0, 0.0, 0.0};

            // float RotationMatrix[3][3];
            // estimatorKalmanGetEstimatedRot((float *)RotationMatrix);
            // float tag_or_versor_calibrated[3] = {RotationMatrix[0][2], RotationMatrix[1][2], RotationMatrix[2][2]};
            float tag_or_versor_calibrated[3] = {0.0, 0.0, 1.0};

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
            // PredictedVoltages[0] = V_rx_1;
            // PredictedVoltages[1] = V_rx_2;
            // PredictedVoltages[2] = V_rx_3;
            // PredictedVoltages[3] = V_rx_4;

            CG_a1 = meanData_a1 / V_rx_1;
            CG_a2 = meanData_a2 / V_rx_2;
            CG_a3 = meanData_a3 / V_rx_3;
            CG_a4 = meanData_a4 / V_rx_4;
            DEBUG_PRINT("CG_a1 = %f\n", (double)CG_a1);
            DEBUG_PRINT("CG_a2 = %f\n", (double)CG_a2);
            DEBUG_PRINT("CG_a3 = %f\n", (double)CG_a3);
            DEBUG_PRINT("CG_a4 = %f\n", (double)CG_a4);
        }
        else
        {
            // ------------------------------------ RUNNING CASE ------------------------------------
            // start_timer();     // start the timer.
            // it1 = get_timer(); // store current cycle-count in a local

            // reset the kalman filter only once after the calibration
            if (currentCalibrationTick == CALIBRATION_TIC_VALUE + 1)
            {
                DEBUG_PRINT("currentCalibrationTick = %f\n", currentCalibrationTick);

                DEBUG_PRINT("Resetting the Kalman filte after calibrationr\n");
                paramSetInt(paramGetVarId("kalman", "resetEstimation"), 1);
                currentCalibrationTick = currentCalibrationTick + 1;
                // paramSetInt(paramGetVarId("kalman", "resetEstimation"), 0);
            }

            // == == == == == == == == == == == == == == == == == == = ADAPTIVE STD COMPUTING -- NOT USED TO NOT OVERPARAMETRIZE THE PROBLEM == == == == == == == == == == == == == == == == == == == == == == == == == =

            if (UseAdaptiveSTD == 1 && number_of_samples_for_DynamicStd < window_size)
            {
                // setting the default std
                Nero_std = voltAnchor->stdDev[0];
                Giallo_std = voltAnchor->stdDev[1];
                Grigio_std = voltAnchor->stdDev[2];
                Rosso_std = voltAnchor->stdDev[3];

                // accumulating the data for the std deviation
                Nero_std_data[number_of_samples_for_DynamicStd] = voltAnchor->measuredVolt[0];
                Giallo_std_data[number_of_samples_for_DynamicStd] = voltAnchor->measuredVolt[1];
                Grigio_std_data[number_of_samples_for_DynamicStd] = voltAnchor->measuredVolt[2];
                Rosso_std_data[number_of_samples_for_DynamicStd] = voltAnchor->measuredVolt[3];

                number_of_samples_for_DynamicStd += 1;
            }
            if (UseAdaptiveSTD == 1 && number_of_samples_for_DynamicStd == window_size)
            {
                // removing the first value
                for (int i = 0; i < window_size - 1; i++)
                {
                    Nero_std_data[i] = Nero_std_data[i + 1];
                    Giallo_std_data[i] = Giallo_std_data[i + 1];
                    Grigio_std_data[i] = Grigio_std_data[i + 1];
                    Rosso_std_data[i] = Rosso_std_data[i + 1];
                }
                // adding the new value
                Nero_std_data[window_size - 1] = voltAnchor->measuredVolt[0];
                Giallo_std_data[window_size - 1] = voltAnchor->measuredVolt[1];
                Grigio_std_data[window_size - 1] = voltAnchor->measuredVolt[2];
                Rosso_std_data[window_size - 1] = voltAnchor->measuredVolt[3];

                // compute the std deviation
                Nero_std = computeSTD(&Nero_std_data[0], window_size);
                Giallo_std = computeSTD(&Giallo_std_data[0], window_size);
                Grigio_std = computeSTD(&Grigio_std_data[0], window_size);
                Rosso_std = computeSTD(&Rosso_std_data[0], window_size);
            }

            // computing the B field for each of the 4 anchors
            point_t cfPosP;
            estimatorKalmanGetEstimatedPos(&cfPosP);
            float cfPos[3] = {cfPosP.x, cfPosP.y, cfPosP.z};

            float tag_pos_predicted[3] = {cfPos[0], cfPos[1], cfPos[2]};
            // float tag_pos_predicted[3] = {0.0, 0.0, 0.0};

            float RotationMatrix[3][3];
            // it would be the product between the rotation matrix and the initial [0,0,1] versor
            estimatorKalmanGetEstimatedRot((float *)RotationMatrix);
            float tag_or_versor[3] = {RotationMatrix[0][2], RotationMatrix[1][2], RotationMatrix[2][2]};
            // float tag_or_versor[3] = {0.0, 0.0, 1.0};

            float anchor_1_pose[3] = {voltAnchor->x[0], voltAnchor->y[0], voltAnchor->z[0]};
            float anchor_2_pose[3] = {voltAnchor->x[1], voltAnchor->y[1], voltAnchor->z[1]};
            float anchor_3_pose[3] = {voltAnchor->x[2], voltAnchor->y[2], voltAnchor->z[2]};
            float anchor_4_pose[3] = {voltAnchor->x[3], voltAnchor->y[3], voltAnchor->z[3]};

            get_B_field_for_a_Anchor(anchor_1_pose, tag_pos_predicted, tag_or_versor, B_field_vector_1);
            get_B_field_for_a_Anchor(anchor_2_pose, tag_pos_predicted, tag_or_versor, B_field_vector_2);
            get_B_field_for_a_Anchor(anchor_3_pose, tag_pos_predicted, tag_or_versor, B_field_vector_3);
            get_B_field_for_a_Anchor(anchor_4_pose, tag_pos_predicted, tag_or_versor, B_field_vector_4);

            // computing the V_rx for each of the 4 anchors
            float V_rx_1 = V_from_B(B_field_vector_1, tag_or_versor, voltAnchor->resonanceFrequency[0], voltAnchor->GainValue);
            float V_rx_2 = V_from_B(B_field_vector_2, tag_or_versor, voltAnchor->resonanceFrequency[1], voltAnchor->GainValue);
            float V_rx_3 = V_from_B(B_field_vector_3, tag_or_versor, voltAnchor->resonanceFrequency[2], voltAnchor->GainValue);
            float V_rx_4 = V_from_B(B_field_vector_4, tag_or_versor, voltAnchor->resonanceFrequency[3], voltAnchor->GainValue);

            // computing the V_rx_derivate for each of the 4 anchors

            V_rx_derivate_x_function(
                tag_pos_predicted[0], tag_pos_predicted[1], tag_pos_predicted[2],
                tag_or_versor[0], tag_or_versor[1], tag_or_versor[2],
                N_WOUNDS, CURRENT, COIL_SURFACE, voltAnchor->GainValue,
                voltAnchor->x[0], voltAnchor->y[0], voltAnchor->z[0],
                voltAnchor->x[1], voltAnchor->y[1], voltAnchor->z[1],
                voltAnchor->x[2], voltAnchor->y[2], voltAnchor->z[2],
                voltAnchor->x[3], voltAnchor->y[3], voltAnchor->z[3],
                voltAnchor->resonanceFrequency[0], voltAnchor->resonanceFrequency[1], voltAnchor->resonanceFrequency[2], voltAnchor->resonanceFrequency[3],
                V_rx_derivative_x);

            V_rx_derivate_y_function(
                tag_pos_predicted[0], tag_pos_predicted[1], tag_pos_predicted[2],
                tag_or_versor[0], tag_or_versor[1], tag_or_versor[2],
                N_WOUNDS, CURRENT, COIL_SURFACE, voltAnchor->GainValue,
                voltAnchor->x[0], voltAnchor->y[0], voltAnchor->z[0],
                voltAnchor->x[1], voltAnchor->y[1], voltAnchor->z[1],
                voltAnchor->x[2], voltAnchor->y[2], voltAnchor->z[2],
                voltAnchor->x[3], voltAnchor->y[3], voltAnchor->z[3],
                voltAnchor->resonanceFrequency[0], voltAnchor->resonanceFrequency[1], voltAnchor->resonanceFrequency[2], voltAnchor->resonanceFrequency[3],
                V_rx_derivative_y);

            V_rx_derivate_z_function(
                tag_pos_predicted[0], tag_pos_predicted[1], tag_pos_predicted[2],
                tag_or_versor[0], tag_or_versor[1], tag_or_versor[2],
                N_WOUNDS, CURRENT, COIL_SURFACE, voltAnchor->GainValue,
                voltAnchor->x[0], voltAnchor->y[0], voltAnchor->z[0],
                voltAnchor->x[1], voltAnchor->y[1], voltAnchor->z[1],
                voltAnchor->x[2], voltAnchor->y[2], voltAnchor->z[2],
                voltAnchor->x[3], voltAnchor->y[3], voltAnchor->z[3],
                voltAnchor->resonanceFrequency[0], voltAnchor->resonanceFrequency[1], voltAnchor->resonanceFrequency[2], voltAnchor->resonanceFrequency[3],
                V_rx_derivative_z);

            // it2 = get_timer() - it1; // Derive the cycle-count difference
            // stop_timer();

            // create the matrix to update the kalman matrix
            float h_1[KC_STATE_DIM] = {0};
            arm_matrix_instance_f32 H_1 = {1, KC_STATE_DIM, h_1};
            h_1[KC_STATE_X] = V_rx_derivative_x[0];
            h_1[KC_STATE_Y] = V_rx_derivative_y[0];
            h_1[KC_STATE_Z] = V_rx_derivative_z[0];

            float h_2[KC_STATE_DIM] = {0};
            arm_matrix_instance_f32 H_2 = {1, KC_STATE_DIM, h_2};
            h_2[KC_STATE_X] = V_rx_derivative_x[1];
            h_2[KC_STATE_Y] = V_rx_derivative_y[1];
            h_2[KC_STATE_Z] = V_rx_derivative_z[1];

            float h_3[KC_STATE_DIM] = {0};
            arm_matrix_instance_f32 H_3 = {1, KC_STATE_DIM, h_3};
            h_3[KC_STATE_X] = V_rx_derivative_x[2];
            h_3[KC_STATE_Y] = V_rx_derivative_y[2];
            h_3[KC_STATE_Z] = V_rx_derivative_z[2];

            float h_4[KC_STATE_DIM] = {0};
            arm_matrix_instance_f32 H_4 = {1, KC_STATE_DIM, h_4};
            h_4[KC_STATE_X] = V_rx_derivative_x[3];
            h_4[KC_STATE_Y] = V_rx_derivative_y[3];
            h_4[KC_STATE_Z] = V_rx_derivative_z[3];

            float error_anchor1 = (voltAnchor->measuredVolt[0] / CG_a1) - V_rx_1; // / 0.995941558f;
            float error_anchor2 = (voltAnchor->measuredVolt[1] / CG_a2) - V_rx_2;
            float error_anchor3 = (voltAnchor->measuredVolt[2] / CG_a3) - V_rx_3;
            float error_anchor4 = (voltAnchor->measuredVolt[3] / CG_a4) - V_rx_4;

            // ------------------- Logging -------------------
            // Predicted Voltages
            PredictedVoltages[0] = V_rx_1;
            PredictedVoltages[1] = V_rx_2;
            PredictedVoltages[2] = V_rx_3;
            PredictedVoltages[3] = V_rx_4;

            // Measured Voltages
            MeasuredVoltages[0] = (voltAnchor->measuredVolt[0] / CG_a1);
            MeasuredVoltages[1] = (voltAnchor->measuredVolt[1] / CG_a2);
            MeasuredVoltages[2] = (voltAnchor->measuredVolt[2] / CG_a3);
            MeasuredVoltages[3] = (voltAnchor->measuredVolt[3] / CG_a4);

            if (UseAdaptiveSTD == 1)
            {
                //  ------------    ------- Case with adaptive std -------------------
                // DEBUG_PRINT("Adaptive\n");
                kalmanCoreScalarUpdate(this, &H_1, error_anchor1, Nero_std);
                kalmanCoreScalarUpdate(this, &H_2, error_anchor2, Giallo_std);
                kalmanCoreScalarUpdate(this, &H_3, error_anchor3, Grigio_std);
                kalmanCoreScalarUpdate(this, &H_4, error_anchor4, Rosso_std);
            }
            else
            {
                // ------------------- Case without adaptive std -------------------
                // DEBUG_PRINT("NON \n");
                // if (printFlagSTD_Type == 0)
                // {
                //     DEBUG_PRINT("NON Adaptive\n");
                //     printFlagSTD_Type = 1;
                // }

                // kalmanCoreScalarUpdate(this, &H_1, error_anchor1, voltAnchor->stdDev[0]);
                // kalmanCoreScalarUpdate(this, &H_2, error_anchor2, voltAnchor->stdDev[1]);
                // kalmanCoreScalarUpdate(this, &H_3, error_anchor3, voltAnchor->stdDev[2]);
                // kalmanCoreScalarUpdate(this, &H_4, error_anchor4, voltAnchor->stdDev[3]);
            }
        }
    }
}

// LOG_GROUP_START(Dipole_Model)
// // // LOG_ADD(LOG_UINT32, CPUCycle, &it2)

// // Measured Voltages
// LOG_ADD(LOG_FLOAT, M_V1, &MeasuredVoltages[0])
// LOG_ADD(LOG_FLOAT, M_V2, &MeasuredVoltages[1])
// LOG_ADD(LOG_FLOAT, M_V3, &MeasuredVoltages[2])
// LOG_ADD(LOG_FLOAT, M_V4, &MeasuredVoltages[3])

// // Predicted Voltages
// LOG_ADD(LOG_FLOAT, P_V1, &PredictedVoltages[0])
// LOG_ADD(LOG_FLOAT, P_V2, &PredictedVoltages[1])
// LOG_ADD(LOG_FLOAT, P_V3, &PredictedVoltages[2])
// LOG_ADD(LOG_FLOAT, P_V4, &PredictedVoltages[3])

// LOG_GROUP_STOP(Dipole_Model)

PARAM_GROUP_START(Dipole_Params)
PARAM_ADD(PARAM_FLOAT, calibTic, &currentCalibrationTick)
PARAM_GROUP_STOP(Dipole_Params)