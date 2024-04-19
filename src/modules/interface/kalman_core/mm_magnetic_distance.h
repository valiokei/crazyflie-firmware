#pragma once
#include "kalman_core.h"

float V_model(kalmanCoreData_t *this,
              voltMeasurement_t *voltAnchor);
