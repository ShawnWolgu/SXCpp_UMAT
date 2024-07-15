#pragma once

#include <Eigen/Eigen>

using namespace Eigen;
using namespace std;

Matrix<double, 6, 6> Jacobian_Matrix(double* ddsdde);
void ddsdde_from_matrix(Matrix<double, 6, 6> Jacobian, double* ddsdde);
