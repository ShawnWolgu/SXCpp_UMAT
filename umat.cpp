#include <Eigen/Eigen>
#include "include/sxcpp.h"
#include "include/func.h"
#include "include/input.h"
#include "include/grain.h"
#include "include/pmode.h"
#include "include/slip.h"
#include <math.h>
#include <stdlib.h>

using namespace std;
using namespace Eigen;

//global variables
int flag_harden, total_mode_num;
Matrix3d lattice_vec;
Matrix6d elastic_modulus_ref;
Matrix6d strain_modi_tensor{
    {1, 0, 0, 0, 0, 0},
    {0, 1, 0, 0, 0, 0},
    {0, 0, 1, 0, 0, 0},
    {0, 0, 0, 2, 0, 0},
    {0, 0, 0, 0, 2, 0},
    {0, 0, 0, 0, 0, 2},
};
LatentMat lat_hard_mat = LatentMat::Identity();
LatentMat interaction_mat = LatentMat::Identity();
PMode* mode_sys[MAX_MODE_NUM] = {nullptr};
Slip slip_array[MAX_MODE_NUM] = {Slip()};

extern "C" void getoutdir(char* outdir, int* lenoutdir, int len);

extern "C" void umat(double* stress, double* statev, double* ddsdde, double* sse, double* spd,
	double* scd, double* rpl, double* ddsddt, double* drplde, double* drpldt,
	double* stran, double* dstran, double* time, double* dtime, double* temp,
	double* dtemp, double* predef, double* dpred, char* cmname, int* ndi,
	int* nshr, int* ntens, int* nstatv, double* props, int* nprops,
	double* coords, double* drot, double* pnewdt, double* celent, double* dfgrd0,
	double* dfgrd1, int* noel, int* npt, int* layer, int* kspt,
	int* kstep, int* kinc, short cmname_len){
    int umat_state = 0;
    // Start of the umat function
    if (*kstep == 1 && *kinc == 0){
        // Initialize the state variables

        /* char temp[200]; */
        /* int lenoutdir; */
        /* getoutdir(temp, &lenoutdir, 200); */
        /* processPath(temp, "\\param.txt"); */
        /* ifstream cij(temp); */

        ifstream cij("param.txt");
        statev[0] = 0;
        statev[1] = 0;
        statev[2] = 0;
        statev[3] = statev[0]; statev[4] = statev[1]; statev[5] = statev[2]; // save the initial euler angle.
        if (total_mode_num == 0){
            elastic_modulus_ref = read_elastic(cij);
            lattice_vec = read_lattice(cij);
            read_pmodes(cij);
            for (int pmode_id = 0; pmode_id < total_mode_num; pmode_id++){
                mode_sys[pmode_id]->cal_shear_modulus(elastic_modulus_ref);
                mode_sys[pmode_id]->initial_statev(statev);
            }
            initialization_interaction();
        }
        else {
            for (int pmode_id = 0; pmode_id < total_mode_num; pmode_id++){
                mode_sys[pmode_id]->initial_statev(statev);
            }
        }
    }
    // Read the state variables
    Matrix3d orientation = Euler_trans(statev[0], statev[1], statev[2]);
    Matrix<double, 6, 6> elastic_modulus = rotate_6d_stiff_modu(elastic_modulus_ref,orientation.transpose());
    Vector6d strain{{stran[0], stran[1], stran[2], stran[5]*0.5, stran[4]*0.5, stran[3]*0.5}};
    Vector6d dstrain_{{dstran[0], dstran[1], dstran[2], dstran[5]*0.5, dstran[4]*0.5, dstran[3]*0.5}};
    Vector6d stress_{{stress[0], stress[1], stress[2], stress[5], stress[4], stress[3]}};
    Matrix3d spin_tensor{
        {drot[0], drot[1], drot[2]}, 
        {drot[3], drot[4], drot[5]}, 
        {drot[6], drot[7], drot[8]}
    };
    spin_tensor = spin_tensor/ *dtime;
    Matrix3d plastic_strain = tensor_trans_order(Vector6d{{statev[6], statev[7], statev[8], statev[9], statev[10], statev[11]}});
    double equiv_plas_strain = statev[12];
    statev[13] = 0.0;//dtemp

    // Jump out if the strain increment is zero: initial stress step
    if (dstrain_.norm() < 1e-20){
        for (int i = 0; i < 6; i++){
            stress[i] += 0;
        }
        ddsdde_from_matrix(elastic_modulus, ddsdde);
        return;
    }
    // dotSigma - WeSigma + SigmaWe + sigma*trace(De) = C De
    // dotSigma - (W-Wp)Sigma + Sigma(W-Wp) + sigma*trace(D-Dp) = C (D-Dp)
    // dotSigma - sigma*trace(Dp) + C Dp = We Sigma - Sigma We - sigma*trace(D) + C D
    Matrix3d spin_plas = Matrix3d::Zero();
    Matrix3d vel_grad_plas = Matrix3d::Zero();
    Matrix3d spin_elas = spin_tensor - spin_plas;
    Matrix3d strain_rate_plas = Matrix3d::Zero();
    Matrix3d strain_rate = tensor_trans_order(dstrain_) / *dtime;
    Matrix6d Cij_pri = get_C_ij_pri(elastic_modulus, stress_);
    Matrix3d stress_3d = tensor_trans_order(stress_);
    Matrix3d jaumann_term = spin_elas * stress_3d - stress_3d * spin_elas;
    Vector6d strain_rate_term = Cij_pri * strain_modi_tensor * tensor_trans_order(strain_rate);
    Vector6d unchanged_term = tensor_trans_order(jaumann_term) + strain_rate_term;
    /* cout << strain_rate.transpose() << endl; */

    // Calculate the stress increment and ddsdde
    Vector6d stress_incr_rate = elastic_modulus * strain_modi_tensor * tensor_trans_order(strain_rate);
    Matrix6d ddp_by_dsigma = Matrix6d::Zero();
    Matrix3d stress_in_iter; Vector6d dp_term; Vector6d F_obj; Matrix6d dF_obj;
    double step_scale = 1.0;
    double F_norm = 1000.0;
    for (int n_iter = 0; n_iter < 20; n_iter++){
        stress_in_iter = stress_3d + tensor_trans_order(stress_incr_rate) * *dtime;
        vel_grad_plas = get_vel_grad_plas(stress_in_iter, orientation, statev, *temp);
        strain_rate_plas = 0.5 * (vel_grad_plas + vel_grad_plas.transpose());
        ddp_by_dsigma = get_dp_grad(stress_in_iter, orientation, statev, *temp);
        dp_term = Cij_pri * strain_modi_tensor * tensor_trans_order(strain_rate_plas);
        F_obj = stress_incr_rate + dp_term - unchanged_term;
        dF_obj = Matrix6d::Identity() + Cij_pri * strain_modi_tensor * ddp_by_dsigma;
        if (isnan(F_obj.norm())){
            /* cout << "[Warning No.1] SXCpp UMAT Warning: NaN values in the stress increment calculation." << endl; */
            umat_state = 1;
            stress_incr_rate = Vector6d::Zero();
            vel_grad_plas = Matrix3d::Zero();
            strain_rate_plas = Matrix3d::Zero();
            ddp_by_dsigma = Matrix6d::Zero();
            break;
        }
        if (F_obj.norm() < 1e-4) break;
        if (n_iter > 8){
            if (F_obj.norm() > 1e-3 || F_obj.norm() > F_norm) {
                step_scale = max(0.1, step_scale * 0.5);
                if (step_scale != 0.1) n_iter -= 1;
            }
            else step_scale = min(1.0, step_scale * 1.5);
        }
        Vector6d dX = dF_obj.inverse() * F_obj;
        stress_incr_rate -= dX * step_scale;
        /* cout << "--" << stress_incr_rate.transpose() << "," << F_obj.norm() << "," << step_scale << endl; */
        F_norm = F_obj.norm();
    }
    if (F_obj.norm() > 1e-3) {
        /* cout << "[Warning No.2] SXCpp UMAT Error: The stress increment calculation did not converge." << endl; */
        umat_state = 2;
    }
    // Update the state variables
    plastic_strain += strain_rate_plas * *dtime;
    Matrix3d elastic_strain = tensor_trans_order(strain) - plastic_strain;
    spin_plas = 0.5 * (vel_grad_plas - vel_grad_plas.transpose());
    spin_elas = spin_tensor - spin_plas;
    orientation = orientation * Rodrigues(spin_elas).transpose(); 
    Vector3d new_euler = Euler_trans(orientation);
    statev[0] = new_euler(0); statev[1] = new_euler(1); statev[2] = new_euler(2);
    Matrix3d stress_grain = tensor_rot_to_CryCoord(stress_in_iter, orientation);
    for (int pmode_id = 0; pmode_id < total_mode_num; pmode_id++){
        mode_sys[pmode_id]->update_ssd(strain_rate, stress_grain, statev, *dtime, *temp);
    }
    /* for (int pmode_id = 0; pmode_id < total_mode_num; pmode_id++){ */
    /*     mode_sys[pmode_id]->update_rho_hard(statev, *dtime, *temp); */
    /* } */
    for (int pmode_id = 0; pmode_id < total_mode_num; pmode_id++){
        mode_sys[pmode_id]->update_status(orientation, elastic_strain, *dtime, *temp, statev);
    }
    stress_ = tensor_trans_order(stress_in_iter);
    Vector6d sigma_out = change_basis_order(stress_);
    for (int i = 0; i < 6; i++){
        stress[i] = sigma_out(i);
    }
    Matrix6d stiffness_eff = elastic_modulus.inverse() + ddp_by_dsigma;
    Matrix6d modulus_eff = stiffness_eff.inverse();
    Matrix6d C_ijkl = change_basis_order(modulus_eff);
    ddsdde_from_matrix(C_ijkl, ddsdde);
}

int main(){
    // Initialize the variables
    double stress[6] = {0};
    double stress_backup[6] = {0};
    double statev[100] = {0};
    double statev_backup[100] = {0};
    double ddsdde[36] = {0};
    double sse = 0, spd = 0, scd = 0, rpl = 0, ddsddt = 0, drplde = 0, drpldt = 0;
    double stran[6] = {0}, dstran[6] = {0};
    double predef = 0, dpred = 0;
    char cmname[1] = {'\0'};
    int ndi = 0, nshr = 0, ntens = 0, nstatv = 13, nprops = 0;
    double props[4] = {10, 10, 10, 0.5};
    double coords[3] = {0}, drot[9] = {0}, pnewdt = 0, celent = 0, dfgrd0[3] = {0}, dfgrd1[3] = {0};
    int noel = 0, npt = 0, layer = 0, kspt = 0, kstep = 1, kinc = 0;
    short cmname_len = 0;

    // Define the deformation rate
    double deformation_rate = 0.01;  // This is just an example value

    // Define the time step
    double time = 0;
    double dtime = 0.01;  // This is just an example value
    // Define the temperature 
    double temp = 300, dtemp = 0;
    //
    // Time-stepping loop
    for (int kinc = 0; kinc < 1000; ++kinc) {
        // Initialize the strain increment
        dstran[0] = deformation_rate * dtime;
        dstran[1] = -deformation_rate * dtime * 0.4866;
        dstran[2] = -deformation_rate * dtime * 0.4866;
        dstran[3] = 0;
        dstran[4] = 0;
        dstran[5] = 0;

        Eigen::MatrixXd F(5, 1); F = Eigen::MatrixXd::Zero(5, 1);
        Eigen::MatrixXd dF(5, 5); dF = Eigen::MatrixXd::Zero(5, 5);

        for (int n_iter = 0; n_iter < 20; n_iter++){
            // Call the umat function
            for (int i = 0; i < 100; i++){
                statev[i] = statev_backup[i];
            }
            for (int i = 0; i < 6; i++){
                stress[i] = stress_backup[i];
            }
            umat(stress, statev, ddsdde, &sse, &spd, &scd, &rpl, &ddsddt, &drplde, &drpldt,
                 stran, dstran, &time, &dtime, &temp, &dtemp, &predef, &dpred, cmname, &ndi,
                 &nshr, &ntens, &nstatv, props, &nprops, coords, drot, &pnewdt, &celent, dfgrd0,
                 dfgrd1, &noel, &npt, &layer, &kspt, &kstep, &kinc, cmname_len);
            for (int i = 1; i < 6; i++){
                F(i-1) = stress[i];
                for (int j = 1; j < 6; j++){
                    dF(i-1, j-1) = ddsdde[6*i + j];
                }
            }
            if (F.norm() < 0.1) break;
            Eigen::MatrixXd dX(5,1);
            dX = dF.inverse() * F;
            for (int i = 0; i < 5; i++){
                dstran[i+1] -= 1* dX(i);
            }
        }
        for (int i = 0; i < 100; i++){
            statev_backup[i] = statev[i];
        }
        for (int i = 0; i < 6; i++){
            stress_backup[i] = stress[i];
        }
        for (int i = 0; i < 6; i++){
            stran[i] += dstran[i];
        }
        // Print the results for this step
        std::cout << stran[0] << "," << stran[1] << "," << stran[2] << "," << stran[3] << "," << stran[4] << "," << stran[5] << ",";
        std::cout << stress[0] << "," << stress[1] << "," << stress[2] << "," << stress[3] << "," << stress[4] << "," << stress[5] << ",";
        std::cout << statev[0] << "," << statev[1] << "," << statev[2] << "," ;
        std::cout << statev[sdv_ind(0,"DD")] << "," << statev[sdv_ind(1,"DD")] << "," << statev[sdv_ind(2,"DD")] << "," << statev[sdv_ind(3,"DD")] << ",";
        std::cout << statev[sdv_ind(4,"DD")] << "," << statev[sdv_ind(5,"DD")] << "," << statev[sdv_ind(6,"DD")] << "," << statev[sdv_ind(7,"DD")] << ",";
        std::cout << statev[sdv_ind(8,"DD")] << "," << statev[sdv_ind(9,"DD")] << "," << statev[sdv_ind(10,"DD")] << "," << statev[sdv_ind(11,"DD")] << endl;
    }
    // Explicitly call destructors for cleanup
    return 0;
}
