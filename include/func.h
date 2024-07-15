#pragma once

#include <Eigen/Eigen>
#include <iostream>
#include <fstream>

using namespace Eigen;
using namespace std;

using Eigen::Matrix3d, Eigen::Vector3d, Eigen::Matrix, Eigen::MatrixXd, Eigen::VectorXd, Eigen::all, Eigen::last;
typedef Matrix<double, 6, 1> Vector6d;
typedef Matrix<double, 6, 6> Matrix6d;

Matrix<double, 6, 6> Jacobian_Matrix(double* ddsdde);
void ddsdde_from_matrix(Matrix<double, 6, 6> Jacobian, double* ddsdde);
Matrix<double, 6, 6> calc_elas_modulus(double* props);
Matrix<double, 6, 6> read_elastic(ifstream &is);
Matrix3d Euler_trans(double phi1, double Phi, double phi2);
Matrix6d rotate_6d_stiff_modu(Matrix6d modulus, Matrix3d rotate_matrix);
Matrix6d cal_rotation_trans_6d_for_stiff(Matrix3d M);
