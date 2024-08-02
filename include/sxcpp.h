#pragma once

#define MAX_MODE_NUM 20
#define MAX_HARDEN_NUM 20
#define MAX_LATENT_NUM 8

#include <Eigen/Eigen>
#include <iostream>
#include <fstream>
#include <cstring>
#include <memory>
#include <array>
#include <new>

using namespace Eigen;
using namespace std;

#define pi 3.14159265358979323846
#define k_boltzmann 1.380649e-23
#define eV_to_J 1.60217662e-19
#define MPa_to_Pa 1e6

typedef Matrix<double, 6, 1> Vector6d;
typedef Matrix<double, 6, 6> Matrix6d;
typedef Matrix<double, MAX_HARDEN_NUM, 1> HardenVec;
typedef Matrix<double, MAX_LATENT_NUM, 1> LatentVec;
typedef Matrix<double, MAX_MODE_NUM, MAX_MODE_NUM> LatentMat;

class Grain;
class PMode;
class Slip;
class Twin;
enum mode_type {slip, twin, undefined};

extern int flag_harden, total_mode_num;
extern Matrix3d lattice_vec;
extern Matrix6d elastic_modulus_ref;
extern Matrix6d strain_modi_tensor;
extern LatentMat lat_hard_mat;
extern LatentMat interaction_mat;
extern array<PMode*, MAX_MODE_NUM> mode_sys;
extern char slip_memory[];
extern Slip* slip_pool[];
