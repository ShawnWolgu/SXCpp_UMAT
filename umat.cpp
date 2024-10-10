#include <Eigen/Eigen>
#include "include/sxcpp.h"
#include "include/func.h"
#include "include/input.h"
#include "include/grain.h"
#include "include/pmode.h"
#include "include/slip.h"
#include <math.h>
#include <stdlib.h>
#include <chrono>
#include <thread>

using namespace std;
using namespace Eigen;

//global variables
int flag_harden, total_mode_num;
bool lock_read=false, lock_wait=true;
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
    double temperature = *temp;
    if (temperature == 0) temperature = 300;
    // Start of the umat function
    if (*kstep == 1 && *kinc < 2){
        // Initialize the state variables
        char temp[200];
        int lenoutdir;
        getoutdir(temp, &lenoutdir, 200);
        processPath(temp, "\\param.txt");
        const char* file_name = temp;

        /* const char* file_name = "param.txt"; */
        /* ifstream cij("param.txt"); */
        statev[0] = 0;
        statev[1] = 0;
        statev[2] = 0;
        statev[3] = statev[0]; statev[4] = statev[1]; statev[5] = statev[2]; // save the initial euler angle.
        if (total_mode_num == 0 && lock_read == false){
            lock_read = true;
            elastic_modulus_ref = read_elastic(file_name);
            lattice_vec = read_lattice(file_name);
            read_pmodes(file_name);
            for (int pmode_id = 0; pmode_id < total_mode_num; pmode_id++){
                mode_sys[pmode_id]->cal_shear_modulus(elastic_modulus_ref);
                mode_sys[pmode_id]->initial_statev(statev);
            }
            initialization_interaction();
            lock_wait = false;
        }
        else {
            if (lock_wait){
                this_thread::sleep_for(chrono::milliseconds(10));
            }
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
    Matrix3d plastic_strain = tensor_trans_order(Vector6d{{statev[6], statev[7], statev[8], statev[9], statev[10], statev[11]}});
    double equiv_plas_strain = statev[12];
    statev[13] = 0.0;//dtemp
    // Jump out if the strain increment is zero: initial stress step
    if (dstrain_.norm() < 1e-20){
        if (*kstep == 1){
            for (int i = 0; i < 6; i++){
                stress[i] += 0;
            }
            ddsdde_from_matrix(elastic_modulus, ddsdde);
            return;
        }
        else{
            dstrain_[0] = 1e-20;
            dstrain_[1] = 1e-20;
            dstrain_[2] = 1e-20;
        }
    }
    // dotSigma - WeSigma + SigmaWe + sigma*trace(De) = C De
    // dotSigma - (W-Wp)Sigma + Sigma(W-Wp) + sigma*trace(D-Dp) = C (D-Dp)
    // dotSigma = C'D + S'W - C'Dp - S'Wp
    Matrix6d C_pri = get_C_ij_pri(elastic_modulus, stress_); //C'
    Matrix6x3d Sigma_Jaumann = get_Sigma_jaumann(stress_); //S'
    Matrix3d strech = tensor_trans_order(dstrain_) / *dtime;
    Matrix3d spin = drot_to_spin(drot, dtime);
    Vector6d unchange_term = C_pri * strain_modi_tensor * tensor_trans_order(strech) +
                             Sigma_Jaumann * tensor_trans_order_spin(spin);

    Matrix3d spin_plas = Matrix3d::Zero();
    Matrix3d vel_grad_plas = Matrix3d::Zero();
    Matrix3d spin_elas = spin - spin_plas;
    Matrix3d strech_plas = Matrix3d::Zero();

    //iteration: Newton-Raphson method
    int iteration_num = 0;
    double step_scale = 1.0, F_norm = 1000.0;
    Vector6d stress_incr_rate = elastic_modulus * strain_modi_tensor * tensor_trans_order(strech);
    Matrix3d stress_3d = tensor_trans_order(stress_);
    Matrix6d ddp_by_dsigma = Matrix6d::Zero();
    Matrix3x6d dwp_by_dsigma = Matrix3x6d::Zero();
    Matrix3d stress_in_iter; Vector6d dp_term, wp_term, F_obj, dX; Matrix6d dF_obj;
    do{
        for (int n_iter = 0; n_iter < 5; n_iter++) {
            stress_in_iter = stress_3d + tensor_trans_order(stress_incr_rate) * *dtime;
            vel_grad_plas = get_vel_grad_plas(stress_in_iter, orientation, statev, temperature);
            strech_plas = 0.5 * (vel_grad_plas + vel_grad_plas.transpose());
            spin_plas = 0.5 * (vel_grad_plas - vel_grad_plas.transpose());
            dp_term = C_pri * strain_modi_tensor * tensor_trans_order(strech_plas);
            wp_term = Sigma_Jaumann * tensor_trans_order_spin(spin_plas);
            F_obj = stress_incr_rate + dp_term + wp_term - unchange_term;

            ddp_by_dsigma = get_dp_grad(stress_in_iter, orientation, statev, temperature);
            dwp_by_dsigma = get_wp_grad(orientation, statev);
            dF_obj = Matrix6d::Identity() + C_pri * strain_modi_tensor * ddp_by_dsigma + Sigma_Jaumann * dwp_by_dsigma;

            F_norm = F_obj.norm();
            if (std::isnan(F_norm)) {
                if (*noel == 2000 && *npt == 1){
                    cout << "NAN detected in the iteration process." << endl;
                    cout << "KSTEP: " << *kstep << " KINC: " << *kinc << endl;
                    cout << "Strecth: " << strech << endl;
                    cout << "Spin: " << spin << endl;
                    cout << "Stress: " << stress_3d << endl;
                    cout << "Stress_incr_rate: " << stress_incr_rate << endl;
                    for (int pmode_id = 0; pmode_id < total_mode_num; pmode_id++){
                        int statev_id = sdv_ind(pmode_id, "SSR");
                        cout << "Mode " << pmode_id << " SSR: " << statev[statev_id] << endl;
                    }
                }
                umat_state = 1;
                stress_incr_rate = elastic_modulus * strain_modi_tensor * tensor_trans_order(strech);
                vel_grad_plas = Matrix3d::Zero();
                strech_plas = Matrix3d::Zero();
                ddp_by_dsigma = Matrix6d::Zero();
                break;
            }
            if (F_norm < CRITERION_CONV) break;
            dX = dF_obj.inverse() * F_obj;
            stress_incr_rate -= dX;
            if (dX.norm() / stress_incr_rate.norm() < 1e-4) break;
            /* cout << "NR iter: " << n_iter << " F_norm: " << F_norm << endl; */
        }
        if (F_norm < CRITERION_CONV) break;
        if (umat_state != 0) break;

        bool downhill_converged = false;
        double initial_F_norm = F_obj.norm();
        double current_F_norm = initial_F_norm;
        step_scale = 1;
        for (int downhill_iter = 0; downhill_iter < 20; downhill_iter++) {
            if (current_F_norm < 1e-4) {
                downhill_converged = true;
                break;
            }
            stress_in_iter = stress_3d + tensor_trans_order(stress_incr_rate) * *dtime;
            vel_grad_plas = get_vel_grad_plas(stress_in_iter, orientation, statev, temperature);
            strech_plas = 0.5 * (vel_grad_plas + vel_grad_plas.transpose());
            spin_plas = 0.5 * (vel_grad_plas - vel_grad_plas.transpose());
            dp_term = C_pri * strain_modi_tensor * tensor_trans_order(strech_plas);
            wp_term = Sigma_Jaumann * tensor_trans_order_spin(spin_plas);
            F_obj = stress_incr_rate + dp_term + wp_term - unchange_term;

            ddp_by_dsigma = get_dp_grad(stress_in_iter, orientation, statev, temperature);
            dwp_by_dsigma = get_wp_grad(orientation, statev);
            dF_obj = Matrix6d::Identity() + C_pri * strain_modi_tensor * ddp_by_dsigma + Sigma_Jaumann * dwp_by_dsigma;

            current_F_norm = F_obj.norm();
            if (std::isnan(current_F_norm)) {
                if (*noel == 2000 && *npt == 1){
                    cout << "NAN detected in the iteration process." << endl;
                    cout << "KSTEP: " << *kstep << " KINC: " << *kinc << endl;
                    cout << "Strecth: " << strech << endl;
                    cout << "Spin: " << spin << endl;
                    cout << "Stress: " << stress_3d << endl;
                    cout << "Stress_incr_rate: " << stress_incr_rate << endl;
                    for (int pmode_id = 0; pmode_id < total_mode_num; pmode_id++){
                        int statev_id = sdv_ind(pmode_id, "SSR");
                        cout << "Mode " << pmode_id << " SSR: " << statev[statev_id] << endl;
                    }
                }
                umat_state = 1;
                stress_incr_rate = elastic_modulus * strain_modi_tensor * tensor_trans_order(strech);
                vel_grad_plas = Matrix3d::Zero();
                strech_plas = Matrix3d::Zero();
                ddp_by_dsigma = Matrix6d::Zero();
                break;
            }
            if (current_F_norm < initial_F_norm) {
                step_scale = std::min(1.0, step_scale * 1.2);
                initial_F_norm = current_F_norm;
                dX = dF_obj.inverse() * F_obj;
                stress_incr_rate -= dX * step_scale;
            } else {
                stress_incr_rate += dX * step_scale * 0.5;
                step_scale = step_scale * 0.5;
            }
            if (dX.norm()*step_scale / stress_incr_rate.norm() < 1e-4) break;
            /* cout << "Downhill iter: " << downhill_iter << " F_norm: " << current_F_norm << endl; */
        }
        /* cout << F_obj.norm() << endl; */
        if (umat_state != 0) break;
        ++iteration_num;
    } while (F_norm > CRITERION_CONV && iteration_num < MAX_ITER_NUM);

    if (F_obj.norm() > 100*CRITERION_CONV) {
        /* cout << "[Warning No.2] SXCpp UMAT Error: The stress increment calculation did not converge." << endl; */
        /* cout << "F_obj: " << F_obj.norm() << endl; */
        umat_state = 2;
    }

    // Update the state variables
    /* cout << "pe_frac: " << calc_equivalent_value(strain_rate_plas)/calc_equivalent_value(strain_rate) << endl; */
    plastic_strain += strech_plas * *dtime;
    Matrix3d elastic_strain = tensor_trans_order(strain) - plastic_strain;
    spin_plas = 0.5 * (vel_grad_plas - vel_grad_plas.transpose());
    spin_elas = spin - spin_plas;
    orientation = orientation * Rodrigues(spin_elas * *time).transpose(); 
    Vector3d new_euler = Euler_trans(orientation);
    statev[0] = new_euler(0); statev[1] = new_euler(1); statev[2] = new_euler(2);
    Matrix3d stress_grain = tensor_rot_to_CryCoord(stress_in_iter, orientation);
    for (int pmode_id = 0; pmode_id < total_mode_num; pmode_id++){
        mode_sys[pmode_id]->update_ssd(strech, stress_grain, statev, *dtime, temperature);
    }
    /* for (int pmode_id = 0; pmode_id < total_mode_num; pmode_id++){ */
    /*     mode_sys[pmode_id]->update_rho_hard(statev, *dtime, *temp); */
    /* } */
    for (int pmode_id = 0; pmode_id < total_mode_num; pmode_id++){
        mode_sys[pmode_id]->update_status(orientation, elastic_strain, *dtime, temperature, statev);
    }
    stress_ = tensor_trans_order(stress_in_iter);
    Vector6d sigma_out = change_basis_order(stress_);
    for (int i = 0; i < 6; i++){
        stress[i] = sigma_out(i);
    }
    Vector6d pe_vec = tensor_trans_order(plastic_strain);
    statev[6] = pe_vec(0); statev[7] = pe_vec(1); statev[8] = pe_vec(2);
    statev[9] = pe_vec(3); statev[10] = pe_vec(4); statev[11] = pe_vec(5);
    equiv_plas_strain = calc_equivalent_value(plastic_strain);
    statev[12] = equiv_plas_strain;
    statev[13] = 0.0;//dtemp
    //
    // Calculate DDSDDE
    Matrix6d d_deps_plas_d_deps = cal_ddEp_ddE(strech_plas, dstrain_, dtime);
    Matrix3x6d d_wpdt_plas_d_deps = cal_dWdt_ddE(spin_plas, dstrain_, dtime);
    Matrix3x6d d_wdt_d_deps = cal_dWdt_ddE(spin, dstrain_, dtime);
    Matrix6d Jacobian = (Matrix6d::Identity() + C_pri * strain_modi_tensor * ddp_by_dsigma * *dtime + Sigma_Jaumann * dwp_by_dsigma * *dtime).inverse() *
                        (C_pri * strain_modi_tensor);
    /* Matrix6d Jacobian = C_pri * strain_modi_tensor - C_pri * strain_modi_tensor * d_deps_plas_d_deps +  */
    /*                     Sigma_Jaumann * d_wdt_d_deps - Sigma_Jaumann * d_wpdt_plas_d_deps; */
    Jacobian = Jacobian * strain_modi_tensor.inverse();
    Matrix6d modulus_eff = change_basis_order(Jacobian) / (1+dstrain_(0)+dstrain_(1)+dstrain_(2));
    for (int i = 0; i < 6; i++){
        for (int j = 0; j < 3; j++){
            modulus_eff(i,j) = modulus_eff(i,j) - stress_(i);
        }
    }
    ddsdde_from_matrix(modulus_eff, ddsdde);
    return;
}
