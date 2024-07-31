#include "sxcpp.h"
#include "pmode.h"
#include "func.h"

inline Matrix3d get_vel_grad_plas(Matrix3d &stress_3d, Matrix3d &orientation, double* statev, double dtime, double temperature){
    Matrix3d vel_grad_plas = Matrix3d::Zero();
    Matrix3d stress_grain = tensor_rot_to_CryCoord(stress_3d, orientation);
    double shear_strain_rate = 0;
    for (int crt_mode_idx = 0; crt_mode_idx < total_mode_num; ++crt_mode_idx){
        shear_strain_rate = mode_sys[crt_mode_idx]->cal_strain(stress_grain, temperature, statev);
        vel_grad_plas += mode_sys[crt_mode_idx]->dL_tensor(shear_strain_rate, dtime);
        statev[sdv_ind(crt_mode_idx, "SSR")] = shear_strain_rate;
    }
    return tensor_rot_to_RefCoord(vel_grad_plas, orientation);
}

inline Matrix6d get_dp_grad(Matrix3d &stress_3d, Matrix3d &orientation, double* statev, double dtime, double temperature){
    Matrix3d stress_grain = tensor_rot_to_CryCoord(stress_3d, orientation);
    double slip_ddgamma_dtau = 0;
    Matrix6d dp_grad = Matrix6d::Zero();
    for (int crt_mode_idx = 0; crt_mode_idx < total_mode_num; ++crt_mode_idx){
        slip_ddgamma_dtau = mode_sys[crt_mode_idx]->cal_ddgamma_dtau(stress_grain, temperature, statev);
        dp_grad += mode_sys[crt_mode_idx]->ddp_dsigma(slip_ddgamma_dtau, dtime);
        statev[sdv_ind(crt_mode_idx, "slope")] = slip_ddgamma_dtau;
    }
    return rotate_6d_compl_modu(dp_grad, orientation.transpose());
}

inline Matrix6d get_wp_grad(Matrix3d &orientation, double* statev, double dtime){
    double slip_ddgamma_dtau = 0;
    Matrix6d dw_grad = Matrix6d::Zero();
    for (int crt_mode_idx = 0; crt_mode_idx < total_mode_num; ++crt_mode_idx){
        slip_ddgamma_dtau = statev[sdv_ind(crt_mode_idx, "slope")];
        dw_grad += mode_sys[crt_mode_idx]->ddp_dsigma(slip_ddgamma_dtau, dtime);
    }
    return rotate_6d_compl_modu(dw_grad, orientation.transpose());
}
