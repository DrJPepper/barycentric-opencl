#pragma once

#include "util.h"

void initialize(MatrixXd*, string);
MatrixXd p(MatrixXd*, float, float);
MatrixXd Q(MatrixXd *points, float u);
float B(float, int);
Vector3d dQdu(MatrixXd*, float, float);
Vector3d dQdv(MatrixXd*, float, float);
MatrixXd normal(MatrixXd*, float, float);
int kChoosei(int n, int k);
int factorial(int);
