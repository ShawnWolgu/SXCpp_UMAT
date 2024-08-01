#include "sxcpp.h"
#include "pmode.h"
#include "func.h"

inline Matrix3d get_vel_grad_plas(Matrix3d &stress_3d, Matrix3d &orientation, double* statev, double temperature){
    Matrix3d vel_grad_plas = Matrix3d::Zero();
    Matrix3d stress_grain = tensor_rot_to_CryCoord(stress_3d, orientation);
    double shear_strain_rate = 0;
    for (int crt_mode_idx = 0; crt_mode_idx < total_mode_num; ++crt_mode_idx){
        shear_strain_rate = mode_sys[crt_mode_idx]->cal_strain(stress_grain, temperature, statev);
        vel_grad_plas += mode_sys[crt_mode_idx]->dL_tensor(shear_strain_rate);
        statev[sdv_ind(crt_mode_idx, "SSR")] = shear_strain_rate;
    }
    return tensor_rot_to_RefCoord(vel_grad_plas, orientation);
}

inline Matrix6d get_dp_grad(Matrix3d &stress_3d, Matrix3d &orientation, double* statev, double temperature){
    Matrix3d stress_grain = tensor_rot_to_CryCoord(stress_3d, orientation);
    double slip_ddgamma_dtau = 0;
    Matrix6d dp_grad = Matrix6d::Zero();
    for (int crt_mode_idx = 0; crt_mode_idx < total_mode_num; ++crt_mode_idx){
        slip_ddgamma_dtau = mode_sys[crt_mode_idx]->cal_ddgamma_dtau(stress_grain, temperature, statev);
        dp_grad += mode_sys[crt_mode_idx]->ddp_dsigma(slip_ddgamma_dtau);
        statev[sdv_ind(crt_mode_idx, "slope")] = slip_ddgamma_dtau;
    }
    return rotate_6d_compl_modu(dp_grad, orientation.transpose());
}

inline Matrix6d get_wp_grad(Matrix3d &orientation, double* statev){
    double slip_ddgamma_dtau = 0;
    Matrix6d dw_grad = Matrix6d::Zero();
    for (int crt_mode_idx = 0; crt_mode_idx < total_mode_num; ++crt_mode_idx){
        slip_ddgamma_dtau = statev[sdv_ind(crt_mode_idx, "slope")];
        dw_grad += mode_sys[crt_mode_idx]->ddp_dsigma(slip_ddgamma_dtau);
    }
    return rotate_6d_compl_modu(dw_grad, orientation.transpose());
}

inline void initialization_interaction(){
    for (int i = 0; i < total_mode_num; i++){
        PMode& imode = *mode_sys[i];
        for (int j = 0; j < total_mode_num; j++){
            if (i == j) {
                lat_hard_mat(i,j) = 1;
                interaction_mat(i,j) = -1;
            }
            else{
                int mode;
                PMode& jmode = *mode_sys[j];
                if (imode.type == slip && jmode.type == slip) {
                    mode = get_interaction_mode(imode.burgers_vec, imode.plane_norm, jmode.burgers_vec, jmode.plane_norm);
                }
                else if (imode.type == slip && jmode.type == twin) mode = 6;
                else if (imode.type == twin && jmode.type == slip) mode = 0; // twin-slip interaction;
                else if (imode.type == twin && jmode.type == twin) mode = 1; // twin-twin interaction;
                else mode = 0;
                lat_hard_mat(i,j) = imode.latent_params[mode];
                interaction_mat(i,j) = mode;
            }
        }
    }
    cout << "Initialization of interaction matrix is done." << endl;
    cout << "interaction matrix:" << endl;
    cout << interaction_mat << endl;
    cout << "lat_hard_mat:" << endl;
    cout << lat_hard_mat << endl;
}

