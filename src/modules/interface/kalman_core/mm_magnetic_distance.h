#pragma once
#include "kalman_core.h"

void kalmanCoreUpdateWithVolt(kalmanCoreData_t *this,
                              voltMeasurement_t *voltAnchor);
