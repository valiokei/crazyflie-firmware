#pragma once
#include "kalman_core.h"

float kalmanCoreUpdateWithVolt(kalmanCoreData_t *this,
                               voltMeasurement_t *voltAnchor);
