#pragma once
#include "sxcpp.h"
#include "func.h"

class PMode {
public:
    Vector3d burgers_vec, plane_norm;
    int num = -1;
    mode_type type = undefined;
    bool flag_active;
    bool if_defined = false;
    double shear_modulus, rho_init;
    /*
     * [velocity parameters] 
     *  1. MFP control coeffient, 2. reference frequency, 3. activation energy, 4. slip resistance, 5. energy exponent
     *  6. saturated speed, 7. drag coefficient
     *
     * [hardening parameters] 
     *  8. forest hardening coefficient
     *
     * [DD evolution parameters] 
     *  0. SSD_density, 9. nucleation coefficient, 10. multiplication coefficient, 11. drag stress D, 12. reference strain rate, 13. c/g 
     *
     * update_params: 0: burgers, 1: mean_free_path, 2: disl_density_resist, 3: forest_stress
     * 
     * [Voce parameters] 
     * 0. tau_0, 1. tau_1, 2. h_0, 3. h_1, 4. ref_strain_rate, 5. SRS
     *
     * [Twin parameters] 
     * 0. tau_0, 1. tau_1, 2. h_0, 3. h_1, 4. twin_strain, 5. SRS, 6. low_threshold, 7. high_threshold*/
    HardenVec harden_params = HardenVec::Zero();
    LatentVec latent_params = LatentVec::Zero();
    Matrix3d schmidt;
    PMode();
    PMode(int slip_num, Vector6d &slip_info, double *hardens, double *latents, Matrix3d lattice_vec, double f_active);
    void cal_shear_modulus(Matrix6d elastic_modulus);
    double cal_rss(Matrix3d stress_tensor);
    virtual double cal_strain(Matrix3d stress_tensor, double temperature, double* statev) {return 0.0;};
    virtual double cal_ddgamma_dtau(Matrix3d stress_tensor, double temperature, double* statev) {return 0.0;};
    virtual void update_status(Matrix3d orientation, Matrix3d strain_elas, double dtime, double temperature, double *statev) {};
    virtual void initial_statev(double *statev) {};
    /* virtual void update_ssd(Matrix3d dstrain, Matrix3d orientation) {}; */
    Matrix3d dL_tensor(double shear_rate, double dtime);
    Matrix6d ddp_dsigma(double ddgamma_dtau, double dtime);
    Matrix6d dwp_dsigma(double ddgamma_dtau, double dtime);
};

inline PMode::PMode() = default;

inline Matrix3d PMode::dL_tensor(double shear_rate, double dtime) {return schmidt * shear_rate;}

inline Matrix6d PMode::ddp_dsigma(double ddgamma_dtau, double dtime) {
    Vector6d symSchmidt_6d = strain_modi_tensor * tensor_trans_order((Matrix3d)(0.5*(schmidt+schmidt.transpose())));
    Matrix6d symSchmidt_66 = symSchmidt_6d * symSchmidt_6d.transpose();
    return symSchmidt_66 * ddgamma_dtau;
}

inline Matrix6d PMode::dwp_dsigma(double ddgamma_dtau, double dtime) {
    Vector6d symSchmidt_6d = strain_modi_tensor * tensor_trans_order((Matrix3d)(0.5*(schmidt+schmidt.transpose())));
    Vector6d asymSchmidt_6d = strain_modi_tensor * tensor_trans_order((Matrix3d)(0.5*(schmidt-schmidt.transpose())));
    Matrix6d symSchmidt_66 = asymSchmidt_6d * symSchmidt_6d.transpose();
    return symSchmidt_66 * ddgamma_dtau;
}

inline double PMode::cal_rss(Matrix3d stress_tensor){
    return (stress_tensor.cwiseProduct(schmidt)).sum();
}

inline void PMode::cal_shear_modulus(Matrix6d elastic_modulus){
    Matrix3d slip_rotation;
    Vector3d trav_direc = burgers_vec.cross(plane_norm);
    slip_rotation << (burgers_vec/burgers_vec.norm()), plane_norm, trav_direc / trav_direc.norm();
    shear_modulus = rotate_6d_stiff_modu(elastic_modulus, slip_rotation.transpose())(3,3);
    if (type == slip){
        cout << "Shear modulus of slip system " << num << " is " << shear_modulus << endl;
    }
    else if (type == twin){
        cout << "Shear modulus of twin system " << num << " is " << shear_modulus << endl;
    }
    else{
        cout << "Shear modulus of deformation system " << num << " is " << shear_modulus << endl;
    }
}

