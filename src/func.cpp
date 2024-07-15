#include "func.h"

Matrix<double, 6, 6> Jacobian_Matrix(double* ddsdde){
    Matrix<double, 6, 6> C;
    for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < 6; j++)
        {
            C(i, j) = ddsdde[i * 6 + j];
        }
    }
    return C;
}

void ddsdde_from_matrix(Matrix<double, 6, 6> Jacobian, double* ddsdde){
    for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < 6; j++)
        {
            ddsdde[i * 6 + j] = Jacobian(i, j);
        }
    }
}
