#include <Eigen/Eigen>
#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include "include/func.h"
#include "include/global.h"

using namespace std;
using namespace Eigen;

//global variables
Matrix<double, 6, 6> elastic_modulus_ref;
extern "C" void getoutdir(char* outdir, int* lenoutdir, int len);

extern "C" void umat(double* stress, double* statev, double* ddsdde, double* sse, double* spd,
	double* scd, double* rpl, double* ddsddt, double* drplde, double* drpldt,
	double* stran, double* dstran, double* time, double* dtime, double* temp,
	double* dtemp, double* predef, double* dpred, char* cmname, int* ndi,
	int* nshr, int* ntens, int* nstatv, double* props, int* nprops,
	double* coords, double* drot, double* pnewdt, double* celent, double* dfgrd0,
	double* dfgrd1, int* noel, int* npt, int* layer, int* kspt,
	int* kstep, int* kinc, short cmname_len){
    // Start of the umat function
    if (*kstep == 1 && *kinc == 0){
        // Initialize the state variables
        char temp[200];
        int lenoutdir;
        getoutdir(temp, &lenoutdir, 200);
        cout << temp << endl;
        ifstream cij("param.txt");
        elastic_modulus_ref = read_elastic(cij);
        statev[0] = 10;
        statev[1] = 10;
        statev[2] = 10;
        cout << elastic_modulus_ref << endl;
    }
    Matrix3d orientation = Euler_trans(statev[0], statev[1], statev[2]);
    Matrix<double, 6, 6> elastic_modulus = rotate_6d_stiff_modu(elastic_modulus_ref,orientation.transpose());
    MatrixXd dstrain_(6, 1);
    MatrixXd stress_(6, 1); //stress in the vector form
    dstrain_ << dstran[0], dstran[1], dstran[2], dstran[3], dstran[4], dstran[5];
    stress_ << stress[0], stress[1], stress[2], stress[3], stress[4], stress[5];
    MatrixXd dstress_ = elastic_modulus * dstrain_;
    for (int i = 0; i < 6; i++){
        stress[i] += dstress_(i);
    }
    ddsdde_from_matrix(elastic_modulus, ddsdde);
}

/* int main(){ */
/*     // Initialize the variables */
/*     double stress[6] = {0}; */
/*     double stress_backup[6] = {0}; */
/*     double statev[13] = {0}; */
/*     double statev_backup[13] = {0}; */
/*     double ddsdde[36] = {0}; */
/*     double sse = 0, spd = 0, scd = 0, rpl = 0, ddsddt = 0, drplde = 0, drpldt = 0; */
/*     double stran[6] = {0}, dstran[6] = {0}, time = 0, dtime = 0, temp = 0, dtemp = 0; */
/*     double predef = 0, dpred = 0; */
/*     char cmname[1] = {'\0'}; */
/*     int ndi = 0, nshr = 0, ntens = 0, nstatv = 13, nprops = 0; */
/*     double props[4] = {10, 10, 10, 0.5}; */
/*     double coords[3] = {0}, drot[3] = {0}, pnewdt = 0, celent = 0, dfgrd0[3] = {0}, dfgrd1[3] = {0}; */
/*     int noel = 0, npt = 0, layer = 0, kspt = 0, kstep = 1, kinc = 0; */
/*     short cmname_len = 0; */
/**/
/*     // Define the deformation rate */
/*     double deformation_rate = 0.01;  // This is just an example value */
/**/
/*     // Define the time step */
/*     double dt = 0.1;  // This is just an example value */
/*     // */
/*     // Time-stepping loop */
/*     for (int kinc = 1; kinc < 11; ++kinc) { */
/*         // Initialize the strain increment */
/*         dstran[0] = deformation_rate * dt; */
/*         dstran[1] = -deformation_rate * dt * 0.4866; */
/*         dstran[2] = -deformation_rate * dt * 0.4866; */
/*         dstran[3] = 0; */
/*         dstran[4] = 0; */
/*         dstran[5] = 0; */
/**/
/*         Eigen::MatrixXd F(5, 1); F = Eigen::MatrixXd::Zero(5, 1); */
/*         Eigen::MatrixXd dF(5, 5); dF = Eigen::MatrixXd::Zero(5, 5); */
/**/
/*         for (int n_iter = 0; n_iter < 10; n_iter++){ */
/*             // Call the umat function */
/*             for (int i = 0; i < 13; i++){ */
/*                 statev[i] = statev_backup[i]; */
/*             } */
/*             for (int i = 0; i < 6; i++){ */
/*                 stress[i] = stress_backup[i]; */
/*             } */
/*             umat(stress, statev, ddsdde, &sse, &spd, &scd, &rpl, &ddsddt, &drplde, &drpldt, */
/*                  stran, dstran, &time, &dtime, &temp, &dtemp, &predef, &dpred, cmname, &ndi, */
/*                  &nshr, &ntens, &nstatv, props, &nprops, coords, drot, &pnewdt, &celent, dfgrd0, */
/*                  dfgrd1, &noel, &npt, &layer, &kspt, &kinc, &kinc, cmname_len); */
/*             for (int i = 1; i < 6; i++){ */
/*                 F(i-1) = stress[i]; */
/*                 for (int j = 1; j < 6; j++){ */
/*                     dF(i-1, j-1) = ddsdde[6*i + j]; */
/*                 } */
/*             } */
/*             if (F.norm() < 0.1) break; */
/*             Eigen::MatrixXd dX(5,1); */
/*             dX = dF.inverse() * F; */
/*             for (int i = 0; i < 5; i++){ */
/*                 dstran[i+1] -= 1* dX(i); */
/*             } */
/*         } */
/*         for (int i = 0; i < 13; i++){ */
/*             statev_backup[i] = statev[i]; */
/*         } */
/*         for (int i = 0; i < 6; i++){ */
/*             stress_backup[i] = stress[i]; */
/*         } */
/*         for (int i = 0; i < 6; i++){ */
/*             stran[i] += dstran[i]; */
/*         } */
/*         // Print the results for this step */
/*         std::cout << stran[0] << "," << stran[1] << "," << stran[2] << "," << stran[3] << "," << stran[4] << "," << stran[5] << ","; */
/*         std::cout << stress[0] << "," << stress[1] << "," << stress[2] << "," << stress[3] << "," << stress[4] << "," << stress[5] << ","; */
/*         std::cout << statev[0] << "," << statev[1] << "," << statev[2] << "," << statev[3] << "," << statev[4] << "," << statev[5] << std::endl; */
/*     } */
/**/
/*     return 0; */
/* } */

