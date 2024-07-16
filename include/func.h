#pragma once

#include <Eigen/Eigen>
#include <iostream>
#include <fstream>
#include <cstring>

using namespace Eigen;
using namespace std;
typedef Matrix<double, 6, 1> Vector6d;
typedef Matrix<double, 6, 6> Matrix6d;

Matrix6d cal_rotation_trans_6d_for_stiff(Matrix3d M);

inline Matrix<double, 6, 6> Jacobian_Matrix(double* ddsdde){
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

inline void ddsdde_from_matrix(Matrix<double, 6, 6> Jacobian, double* ddsdde){
    for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < 6; j++)
        {
            ddsdde[i * 6 + j] = Jacobian(i, j);
        }
    }
}

inline Matrix<double, 6, 6> read_elastic(ifstream &is){
    int row_num = 0, temp_idx = 0;
    double temp[6] = {0,0,0,0,0,0};
    string temp_str;
    Matrix<double, 6, 6> modulus;
    for(;row_num !=6; ++row_num){
        temp_idx = 0;
        getline(is, temp_str);
        cout << temp_str << endl;
        stringstream stream(temp_str);
        while(!stream.eof() && temp_idx!=6) stream >> temp[temp_idx++];
        modulus.row(row_num) << temp[0], temp[1], temp[2], temp[3], temp[4], temp[5];
    }
    /* cout << "Read Elastic Modulus:" << endl << modulus << endl; */
    return modulus;
}

inline Matrix3d Euler_trans(double phi1, double Phi, double phi2){
    Matrix3d Mout;
    double SPH,STH,STM;
    double CPH,CTH,CTM;
    SPH=sin(phi1); CPH=cos(phi1);
    STH=sin(Phi); CTH=cos(Phi);
    STM=sin(phi2); CTM=cos(phi2);
    Mout(0,0)=CTM*CPH-SPH*STM*CTH;
    Mout(1,0)=-STM*CPH-SPH*CTM*CTH;
    Mout(2,0)=SPH*STH;
    Mout(0,1)=CTM*SPH+CPH*STM*CTH;
    Mout(1,1)=-SPH*STM+CPH*CTM*CTH;
    Mout(2,1)=-STH*CPH;
    Mout(0,2)=STH*STM;
    Mout(1,2)=CTM*STH;
    Mout(2,2)=CTH;
    return Mout;
}

inline Matrix6d rotate_6d_stiff_modu(Matrix6d modulus, Matrix3d rotate_matrix){
    Matrix6d M66 = cal_rotation_trans_6d_for_stiff(rotate_matrix);
    return M66*modulus*M66.transpose();
}

inline Matrix6d cal_rotation_trans_6d_for_stiff(Matrix3d M){
   Matrix6d M66;
    double xx,yy,zz,yz,xz,xy,yx,zx,zy;
    xx = M(0,0); 
    yy = M(1,1);
    zz = M(2,2);
    yz = M(1,2);
    xz = M(0,2);
    xy = M(0,1);
	yx = M(1,0);
	zx = M(2,0);
	zy = M(2,1);

    for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
    M66(i,j) = M(i,j)*M(i,j);

    M66(0,3) = 2*xy*xz;
    M66(0,4) = 2*xx*xz;
    M66(0,5) = 2*xx*xy;

    M66(1,3) = 2*yy*yz;
    M66(1,4) = 2*yx*yz;
    M66(1,5) = 2*yx*yy;

    M66(2,3) = 2*zy*zz;
    M66(2,4) = 2*zx*zz;
    M66(2,5) = 2*zx*zy;

    M66(3,0) = yx*zx;
    M66(3,1) = yy*zy;
    M66(3,2) = yz*zz;

    M66(4,0) = xx*zx;
    M66(4,1) = xy*zy;
    M66(4,2) = xz*zz;

    M66(5,0) = xx*yx;
    M66(5,1) = xy*yy;
    M66(5,2) = xz*yz;

    M66(3,3) = yy*zz+yz*zy;
    M66(3,4) = yx*zz+yz*zx;
    M66(3,5) = yx*zy+yy*zx;

    M66(4,3) = xy*zz+xz*zy;
    M66(4,4) = xx*zz+xz*zx;
    M66(4,5) = xx*zy+xy*zx;

    M66(5,3) = xy*yz+xz*yy;
    M66(5,4) = xx*yz+xz*yx;
    M66(5,5) = xx*yy+xy*yx;

	return M66;
}

inline void processPath(char* path, const char* terminationSequence) {
    const size_t terminationLength = strlen(terminationSequence);
    const size_t maxPathLength = 200;
    // Find the length of the valid path part (non-empty part)
    size_t len = 0;
    for (size_t i = 0; i < maxPathLength && path[i] != '\0'; ++i) {
        if (path[i] != ' ') {
            len = i + 1;
        }
    }

    // Add the termination sequence after the valid path
    if (len + terminationLength < 200) { // Ensure there's enough space for the sequence and null terminator
        strcpy(path + len, terminationSequence);
    } else {
        std::cerr << "Error: Not enough space to add termination sequence." << std::endl;
    }
}

