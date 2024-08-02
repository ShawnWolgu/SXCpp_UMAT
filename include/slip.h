#pragma once
#include "func.h"
#include "pmode.h"
#include "sxcpp.h"
#include <array>
#include <memory>

double waiting_time(double stress_eff, double resistance_slip,
                    double act_energy_r, double frequency_r, double energy_expo,
                    double temperature);
double running_time(double stress_eff, double c_drag, double speed_sat,
                    double mean_free_path, double burgers, double temperature);
double waiting_time_grad(double stress_eff, double resistance_slip,
                                 double act_energy_r, double frequency_r,
                                 double energy_expo, double temperature);
double running_time_grad(double stress_eff, double c_drag,
                                 double speed_sat, double mean_free_path,
                                 double burgers, double temperature);

class Slip : public PMode {
public:
    Slip();
    Slip(int slip_num, Vector6d &slip_info, HardenVec &hardens,
         LatentVec &latents, Matrix3d lattice_vec, double f_active);
    double cal_strain(Matrix3d stress_tensor, double temperature, double* statev) override;
    double cal_ddgamma_dtau(Matrix3d stress_tensor, double temperature, double* statev) override;
    void update_status(Matrix3d orientation, Matrix3d strain_elas, double dtime, double temperature, double* statev) override;
    void initial_statev(double* statev) override;
    void update_ssd(Matrix3d strain_rate, Matrix3d stress_grain, double *statev, double dtime, double temperature) override;
    void update_rho_hard(double *statev, double dtime, double temperature) override;
    /* void update_ssd(Matrix3d dstrain, Matrix3d orientation) override; */
    virtual ~Slip() = default;

private:
    double cal_strain_pow(Matrix3d stress_tensor, double* statev);
    double cal_strain_disvel(Matrix3d stress_tensor, double temperature, double* statev);
    double cal_ddgamma_dtau_pow(Matrix3d stress_tensor, double* statev);
    double cal_ddgamma_dtau_disvel(Matrix3d stress_tensor, double temperature, double* statev);
    void update_voce(LatentMat &lat_hard_mat, double dtime, double *statev);
    void update_disvel(LatentMat &lat_hard_mat, double bv,
                       double dtime, double temperature, double* statev);
    double disl_velocity(double rss, double burgers, double temperature, double* statev);
    double disl_velocity_grad(double rss, double burgers, double temperature, double* statev);
};

inline Slip::Slip() = default;

inline Slip::Slip(int slip_num, Vector6d &slip_info, HardenVec &hardens,
                  LatentVec &latents, Matrix3d lattice_vec, double f_active) {
    Vector3d plane_norm_disp;
    num = slip_num, type = slip;
    harden_params = hardens;
    latent_params = latents;
    flag_active = !(f_active < 1e-20);
    for (int temp_idx = 0; temp_idx < 6; ++temp_idx) {
        if (temp_idx < 3)
            plane_norm_disp(temp_idx) = slip_info(temp_idx);
        else {
            if (temp_idx < 6)
                burgers_vec(temp_idx - 3) = slip_info(temp_idx);
        }
    }
    burgers_vec = (burgers_vec.transpose() * lattice_vec).transpose();
    plane_norm = get_plane_norm(plane_norm_disp, lattice_vec);
    schmidt = burgers_vec / burgers_vec.norm() * plane_norm.transpose();
    if (flag_harden == 1) rho_init = harden_params[0] * f_active;
    shear_modulus = 0;
    cout << "Slip system " << num << " is initialized." << endl;
    cout << "Burgers vector of slip system " << num << " is "
        << burgers_vec.transpose() << endl;
    cout << "Normal vector of slip system " << num << " is "
        << plane_norm.transpose() << endl;
    cout << "Harden parameters of slip system " << num << " is ";
    for (int temp_idx = 0; temp_idx < harden_params.size(); ++temp_idx) {
        cout << harden_params[temp_idx] << " ";
    }
    cout << endl;
    cout << "Latent parameters of slip system " << num << " is ";
    for (int temp_idx = 0; temp_idx < latent_params.size(); ++temp_idx) {
        cout << latent_params[temp_idx] << " ";
    }
    cout << endl;
    if_defined = true;
};

inline double Slip::cal_strain(Matrix3d stress_tensor, double temperature, double *statev) {
    /*
   * Select different model for shear strain rate, controlled by flag_harden.
   * 0 : Voce Hardening; 1 : Dislocation Velocity Model;
   * Power model will be used in case 0, while Velocity model used in case 1;
   */
    if (flag_active == false) {
        return 0.0;
    }
    double shear_rate = 0.0;
    switch (flag_harden) {
        case 0:
            shear_rate = cal_strain_pow(stress_tensor, statev);
            break;
        case 1:
            shear_rate = cal_strain_disvel(stress_tensor, temperature, statev);
            break;
        default:
            shear_rate = cal_strain_pow(stress_tensor, statev);
            break;
    }
    return shear_rate;
}

inline double Slip::cal_strain_pow(Matrix3d stress_tensor, double *statev) {
    double shear_rate = 0.0;
    double rss_slip = cal_rss(stress_tensor);
    double crss = statev[sdv_ind(num, "CRSS")];
    double ref_strain_rate = harden_params[4];
    double rate_sen = harden_params[5];
    if (abs(rss_slip) > 0.5 * crss) {
        shear_rate = ref_strain_rate * pow(abs(rss_slip / crss), 1 / rate_sen) * sign(rss_slip);
    }
    return shear_rate;
}

inline double Slip::cal_strain_disvel(Matrix3d stress_tensor,
                                      double temperature, double *statev){
    double shear_rate;
    double burgers = burgers_vec.norm() * 1e-10;
    double rss_slip = cal_rss(stress_tensor);
    double disl_vel = disl_velocity(rss_slip, burgers, temperature, statev);
    double disl_density = statev[sdv_ind(num, "DD")];
    shear_rate = abs(disl_density * burgers * disl_vel) * sign(rss_slip);
    return shear_rate;
}

inline void Slip::update_status(Matrix3d orientation, Matrix3d strain_elas,
                                double dtime, double temperature, double *statev) {
    /* Select different model for shear strain rate, controlled by flag_harden. */
    /* 0 : Voce Hardening; 1 : Dislocation Velocity Model; */
    /* Update Schmidt here. */
    /* Vector3d update_bv = grain.orientation * grain.deform_grad_elas *
   * grain.orient_ref.transpose() * burgers_vec; */
    /* Vector3d update_nv = grain.orientation *
   * grain.deform_grad_elas.inverse().transpose() * grain.orient_ref.transpose()
   * * plane_norm; */
    /* schmidt = update_bv / update_bv.norm() * update_nv.transpose(); */
    Vector3d update_bv = burgers_vec;
    switch (flag_harden) {
        case 0:
            update_voce(lat_hard_mat, dtime, statev);
            break;
        case 1:
            update_disvel(lat_hard_mat, update_bv.norm(), dtime, temperature, statev);
            break;
        default:
            update_voce(lat_hard_mat, dtime, statev);
            break;
    }
}

/* Update crss */
inline void Slip::update_voce(LatentMat &lat_hard_mat, double dtime, double *statev) {
    double Gamma = 0; double crss = statev[sdv_ind(num, "CRSS")];
    for (int crt_mode_idx = 0; crt_mode_idx < total_mode_num; ++crt_mode_idx) {
        Gamma += abs(statev[sdv_ind(crt_mode_idx, "ACC")]);
    }
    double tau_0 = harden_params[0], tau_1 = harden_params[1],
    h_0 = harden_params[2], h_1 = harden_params[3];
    double dtau_by_dGamma =
        h_1 + (abs(h_0 / tau_1) * tau_1 - h_1) * exp(-Gamma * abs(h_0 / tau_1)) +
        abs(h_0 / tau_1) * h_1 * Gamma * exp(-Gamma * abs(h_0 / tau_1));
    for (int crt_mode_idx = 0; crt_mode_idx < total_mode_num; ++crt_mode_idx) {
        crss += abs(statev[sdv_ind(crt_mode_idx, "SSR")]) 
            * dtime * lat_hard_mat(num, mode_sys[crt_mode_idx]->num) * dtau_by_dGamma;
    }
    statev[sdv_ind(num, "CRSS")] = crss;
}

inline void Slip::update_disvel(LatentMat &lat_hard_mat, double bv,
                       double dtime, double temperature, double* statev) {
    /*
     * [velocity parameters] 
     *  1. MFP control coeffient, 2. reference frequency, 3. activation energy, 4. slip resistance, 5. energy exponent
     *  6. saturated speed, 7. drag coefficient
     * [hardening parameters] 
     *  8. forest hardening coefficient
     * [DD evolution parameters] 
     *  0. SSD_density, 9. nucleation coefficient, 10. multiplication coefficient, 11. drag stress D, 12. reference strain rate, 13. c/g 
     * [Twin parameters] 
     * 0. tau_0, 1. tau_1, 2. h_0, 3. h_1, 4. twin_strain, 5. SRS, 6. low_threshold, 7. high_threshold*/
    double burgers = bv * 1e-10, joint_factor = 0.0;
    double c_mfp = harden_params[1], c_forest = harden_params[8], HP_stress = 0;
    double disl_density, disl_density_for, disl_density_resist, joint_density, forest_stress, mean_free_path;
    disl_density = statev[sdv_ind(num, "DD")];
    disl_density_for = disl_density_resist = joint_density = 0;
    for (int crt_mode_idx = 0; crt_mode_idx < total_mode_num; ++crt_mode_idx) {
        disl_density_for += statev[sdv_ind(crt_mode_idx, "DD")];
        disl_density_resist += statev[sdv_ind(crt_mode_idx, "DD")] * lat_hard_mat(num, crt_mode_idx);
        double rho_asterisk = statev[sdv_ind(crt_mode_idx, "DD")] - mode_sys[crt_mode_idx]->rho_init;
        if (crt_mode_idx != num)
            joint_density += lat_hard_mat(num, crt_mode_idx) * sqrt(rho_asterisk) * sqrt(disl_density-rho_init);
    }
    double crss_factor = joint_factor * joint_density + disl_density_resist;
    forest_stress = c_forest * shear_modulus * burgers * sqrt(crss_factor) + HP_stress;
    mean_free_path = c_mfp / sqrt(disl_density_for);
    statev[sdv_ind(num,"mfp")] = mean_free_path; statev[sdv_ind(num,"tauf")] = forest_stress;
}

inline double Slip::cal_ddgamma_dtau(Matrix3d stress_tensor, double temperature, double* statev) {
    /* Select different model for shear strain rate gradient calculation,
   * controlled by flag_harden. */
    /* 0 : Voce Hardening; 1 : Dislocation Velocity Model; */
    /* Update Schmidt here. */
    /* Vector3d update_bv =  burgers_vec; */
    /* Vector3d update_nv =  plane_norm; */
    if (flag_active) {
        statev[sdv_ind(num, "slope")] = 0.0;
        return 0.0;
    }
    switch (flag_harden) {
        case 0:
            statev[sdv_ind(num, "slope")] = cal_ddgamma_dtau_pow(stress_tensor, statev);
            break;
        case 1:
            statev[sdv_ind(num, "slope")] = cal_ddgamma_dtau_disvel(stress_tensor, temperature, statev);
            break;
        default:
            statev[sdv_ind(num, "slope")] = cal_ddgamma_dtau_pow(stress_tensor, statev);
            break;
    }
    return statev[sdv_ind(num, "slope")];
}

inline double Slip::cal_ddgamma_dtau_pow(Matrix3d stress_tensor, double* statev) {
    double ddgamma_dtau = 0.0;
    double rss_slip = cal_rss(stress_tensor);
    double crss = statev[sdv_ind(num, "CRSS")];
    double ref_strain_rate = harden_params[4];
    double rate_sen = harden_params[5];
    if (abs(rss_slip) > 0.5 * crss) {
        ddgamma_dtau = ref_strain_rate *
            pow(abs(rss_slip / crss), 1 / rate_sen - 1) *
            sign(rss_slip) / rate_sen / crss * sign(rss_slip);
    }
    return ddgamma_dtau;
}

inline double Slip::cal_ddgamma_dtau_disvel(Matrix3d stress_tensor, double temperature, double* statev) {
    double burgers = burgers_vec.norm() * 1e-10;
    double rss_slip = cal_rss(stress_tensor);
    double disl_density = statev[sdv_ind(num, "DD")];
    double dvel = disl_velocity_grad(rss_slip, burgers, temperature, statev);
    double ddgamma_dtau = disl_density * burgers * sign(rss_slip) * dvel;
    return ddgamma_dtau;
}

inline void Slip::initial_statev(double* statev){
    if (statev == nullptr) {
        cout << "State variables are not initialized." << endl;
        return;
    }
    if (flag_harden == 0) {
        statev[sdv_ind(num, "SR")] = 0.0;
        statev[sdv_ind(num, "slope")] = 0.0;
        statev[sdv_ind(num, "ACC")] = 0.0;
        statev[sdv_ind(num, "CRSS")] = harden_params[0];
        statev[sdv_ind(num, "custom")] = 0.0;
    } 
    else if (flag_harden == 1) {
        statev[sdv_ind(num, "SR")] = 0.0;
        statev[sdv_ind(num, "slope")] = 0.0;
        statev[sdv_ind(num, "DD")] = rho_init;
        statev[sdv_ind(num, "tauf")] = harden_params[8] * shear_modulus * burgers_vec.norm() * 1e-10 * sqrt(rho_init*total_mode_num);
        statev[sdv_ind(num, "mfp")] = harden_params[1]/sqrt(rho_init*total_mode_num)*total_mode_num;
    }
    else {
        statev[sdv_ind(num, "SR")] = 0.0;
        statev[sdv_ind(num, "slope")] = 0.0;
        statev[sdv_ind(num, "ACC")] = 0.0;
        statev[sdv_ind(num, "CRSS")] = harden_params[0];
        statev[sdv_ind(num, "custom")] = 0.0;
    }
}

inline double Slip::disl_velocity(double rss, double burgers, double temperature, double* statev) {
    /*
   * [velocity parameters]
   *  1. MFP control coeffient, 2. reference frequency, 3. activation energy, 4.
   * slip resistance, 5. energy exponent
   *  6. saturated speed, 7. drag coefficient
   * [hardening parameters]
   *  8. forest hardening coefficient
   * [DD evolution parameters]
   *  0. SSD_density, 9. multiplication coefficient, 10. drag stress D, 11.
   * reference strain rate, 12. c/g
   *
   * update parameters:
   * 0: burgers, 1: mean_free_path, 2: disl_density_resist, 3: forest_stress,
   */
    double frequency_r = harden_params[2], act_energy_r = harden_params[3],
    resistance_slip = harden_params[4], energy_expo = harden_params[5],
    speed_sat = harden_params[6], c_drag = harden_params[7];
    double mean_free_path = statev[sdv_ind(num, "mfp")], t_wait, t_run;
    double forest_stress = statev[sdv_ind(num, "tauf")];
    double stress_eff = abs(rss) - forest_stress;
    if (stress_eff >= 0.0) {
        t_wait = waiting_time(stress_eff, resistance_slip, act_energy_r,
                              frequency_r, energy_expo, temperature);
        t_run = running_time(stress_eff, c_drag, speed_sat, mean_free_path, burgers,
                             temperature);
        return mean_free_path / (t_wait + t_run);
    } else {
        t_wait = waiting_time(stress_eff, resistance_slip, act_energy_r,
                              frequency_r, energy_expo, temperature);
        return mean_free_path / t_wait;
    }
}

inline double Slip::disl_velocity_grad(double rss, double burgers, double temperature, double* statev) {
    /*
   * [velocity parameters]
   *  1. MFP control coeffient, 2. reference frequency, 3. activation energy, 4.
   * slip resistance, 5. energy exponent
   *  6. saturated speed, 7. drag coefficient
   * [hardening parameters]
   *  8. forest hardening coefficient
   * [DD evolution parameters]
   *  0. SSD_density, 9. multiplication coefficient, 10. drag stressD, 11.
   * reference strain rate, 12. c/g
   *
   * update parameters:
   * 0: burgers, 1: mean_free_path, 2: disl_density_resist, 3: forest_stress,
   */
    double frequency_r = harden_params[2], act_energy_r = harden_params[3],
    resistance_slip = harden_params[4], energy_expo = harden_params[5],
    speed_sat = harden_params[6], c_drag = harden_params[7];
    double mean_free_path = statev[sdv_ind(num, "mfp")], t_wait, t_run;
    double forest_stress = statev[sdv_ind(num, "tauf")];
    double stress_eff = abs(rss) - forest_stress;
    double dvel_dtau = 0.0;
    if (stress_eff >= 0.0) {
        t_wait = waiting_time(stress_eff, resistance_slip, act_energy_r,
                              frequency_r, energy_expo, temperature);
        t_run = running_time(stress_eff, c_drag, speed_sat, mean_free_path, burgers,
                             temperature);
        double dtwait_drss = waiting_time_grad(stress_eff, resistance_slip, act_energy_r, frequency_r, energy_expo, temperature);
        double dtrun_drss = running_time_grad(stress_eff, c_drag, speed_sat, mean_free_path, burgers, temperature);
        dvel_dtau = -sign(rss) * mean_free_path / pow((t_wait + t_run), 2) * (dtwait_drss + dtrun_drss);
    } else {
        t_wait = waiting_time(stress_eff, resistance_slip, act_energy_r,
                              frequency_r, energy_expo, temperature);
        double dtwait_drss =
            waiting_time_grad(stress_eff, resistance_slip, act_energy_r,
                              frequency_r, energy_expo, temperature);
        dvel_dtau = -sign(rss) * mean_free_path / pow(t_wait, 2) * dtwait_drss;
    }
    return dvel_dtau;
}

inline double waiting_time(double stress_eff, double resistance_slip,
                           double act_energy_r, double frequency_r,
                           double energy_expo, double temperature) {
    stress_eff = stress_eff * MPa_to_Pa,
        resistance_slip = resistance_slip * MPa_to_Pa,
        act_energy_r = act_energy_r * eV_to_J;
    double act_energy = 0.;
    if (stress_eff >= 0.0)
        act_energy = act_energy_r * (1 - pow((stress_eff / resistance_slip),
                                             energy_expo)); // activation energy
    else
        act_energy = act_energy_r * (1 + pow((abs(stress_eff) / resistance_slip),
                                             energy_expo)); // activation energy
    act_energy =
        min(act_energy, 500 * k_boltzmann *
            temperature); // avoid too large activation energy;
    act_energy =
        max(act_energy, -500 * k_boltzmann *
            temperature); // avoid too small activation energy;
    return 1 / frequency_r * exp(act_energy / (k_boltzmann * temperature));
}

inline double waiting_time_grad(double stress_eff,
                                        double resistance_slip,
                                        double act_energy_r, double frequency_r,
                                        double energy_expo,
                                        double temperature) {
    /* Return a vector: 0. dtw/dtau, 1. tw. */
    stress_eff = stress_eff * MPa_to_Pa,
        resistance_slip = resistance_slip * MPa_to_Pa,
        act_energy_r = act_energy_r * eV_to_J;
    double act_energy = 0.;
    if (stress_eff >= 0.0)
        act_energy = act_energy_r * (1 - pow((stress_eff / resistance_slip),
                                             energy_expo)); // activation energy
    else
        act_energy = act_energy_r * (1 + pow((abs(stress_eff) / resistance_slip),
                                             energy_expo)); // activation energy
    act_energy =
        min(act_energy, 500 * k_boltzmann *
            temperature); // avoid too large activation energy;
    act_energy =
        max(act_energy, -500 * k_boltzmann *
            temperature); // avoid too small activation energy;
    double waiting_time =
        1 / frequency_r * exp(act_energy / (k_boltzmann * temperature));
    double grad_term = 0;
    if (act_energy != 500 * k_boltzmann * temperature &&
        act_energy != -500 * k_boltzmann * temperature)
        grad_term = -act_energy_r * energy_expo *
            pow((abs(stress_eff) / resistance_slip), energy_expo - 1) /
            (k_boltzmann * temperature);
    return grad_term * waiting_time * MPa_to_Pa / resistance_slip;
}

inline double running_time(double stress_eff, double c_drag, double speed_sat,
                           double mean_free_path, double burgers,
                           double temperature) {
    if (stress_eff <= 1e-4) return 1e-40;
    stress_eff = stress_eff * MPa_to_Pa;
    double coeff_B =
        (c_drag * k_boltzmann * temperature) / (speed_sat * pow(burgers, 2));
    /* double coeff_B = (c_drag * k_boltzmann * temperature) / (speed_sat *
   * burgers * mean_free_path); */
    double v_norm = 2 * burgers * stress_eff / (coeff_B * speed_sat);
    double velocity = speed_sat * (sqrt(1 + pow(v_norm, -2)) - 1 / v_norm);
    return mean_free_path / max(velocity, 1e-40);
}

inline double running_time_grad(double stress_eff, double c_drag,
                                        double speed_sat, double mean_free_path,
                                        double burgers, double temperature) {
    /* Return a vector: 0. dtr/dtau, 1. tr. */
    stress_eff = stress_eff * MPa_to_Pa;
    double coeff_B =
        (c_drag * k_boltzmann * temperature) / (speed_sat * burgers * burgers);
    /* double coeff_B = (c_drag * k_boltzmann * temperature) / (speed_sat *
   * burgers * mean_free_path); */
    double v_norm = 2 * burgers * stress_eff / (coeff_B * speed_sat);
    double velocity = speed_sat * (sqrt(1 + 1 / pow(v_norm, 2)) - 1 / v_norm);
    double gradient = -1 * (2 * burgers * mean_free_path) /
        (pow(velocity, 2) * coeff_B) * pow(v_norm, -2) *
        (1 - 1 / (v_norm * sqrt(1 + pow(v_norm, -2))));
    return gradient * MPa_to_Pa;
}

inline void Slip::update_ssd(Matrix3d strain_rate, Matrix3d stress, double *statev, double dtime, double temperature){
    /*
     * [velocity parameters] 
     *  1. MFP control coeffient, 2. reference frequency, 3. activation energy, 4. slip resistance, 5. energy exponent
     *  6. saturated speed, 7. drag coefficient
     * [hardening parameters] 
     *  8. forest hardening coefficient
     * [DD evolution parameters] 
     *  0. SSD_density, 9. nucleation coefficient, 10. nucleation threshold stress, 11. multiplication coefficient
     *  12. drag stress D, 13. reference strain rate, 14. c/g 
     *
     * update parameters:
     * 0: burgers, 1: mean_free_path, 2: disl_density_resist, 3: forest_stress
     */
    if (flag_harden == 0){
        statev[sdv_ind(num,"ACC")] += statev[sdv_ind(num,"SSR")] * dtime;
    }
    if (flag_harden == 1){ 
        double c_forest = harden_params[8], c_nuc = harden_params[9], tau_nuc = harden_params[10],\
               c_multi = harden_params[11], c_annih = 0.,\
               D = harden_params[12] * 1e6, ref_srate = harden_params[13], gg = c_forest/harden_params[14],\
               burgers = burgers_vec.norm() * 1e-10, mfp = statev[sdv_ind(num,"mfp")], \
               forest_stress = statev[sdv_ind(num,"tauf")], SSD_density = statev[sdv_ind(num,"DD")], \
               shear_rate = statev[sdv_ind(num,"SSR")], rho_sat = 0., rho_mov = 0.;
        /* double equi_strain_rate = strain_rate(2,2); */
        /* double equi_strain_rate = calc_first_principal(strain_rate); */
        /* double equi_strain_rate = abs(shear_rate); */
        double equi_strain_rate = calc_equivalent_value(strain_rate);
        rho_sat = c_forest * burgers / gg * (1-k_boltzmann * temperature/D/pow(burgers,3) * log(abs(equi_strain_rate)/ref_srate));
        rho_sat = max(pow(1/rho_sat,2), 0.5*SSD_density);
        if (equi_strain_rate < 1e-20) rho_sat = SSD_density;
        /* double eff_nuc_stress = max(abs(rss) - forest_stress, 0.) - tau_nuc; */
        double eff_nuc_stress = abs(cal_rss(stress)) - tau_nuc;
        double term_nuc = c_nuc * max(eff_nuc_stress,0.) / (shear_modulus * burgers * burgers);
        double term_multi = c_multi / mfp; 
        c_annih = (term_multi + term_nuc) / rho_sat;
        double disloc_incre = (term_multi + term_nuc - c_annih * SSD_density) * abs(shear_rate) * dtime;
        if (disloc_incre > 2 * rho_sat) {
            disloc_incre = 2 * rho_sat - SSD_density; 
        }
        else if (disloc_incre + SSD_density < 0) disloc_incre = -0.1 * SSD_density; 
        SSD_density += disloc_incre;
        rho_mov = SSD_density;
        if(SSD_density < rho_init) rho_init = SSD_density;
        statev[sdv_ind(num,"DD")] = SSD_density;
    }
}

inline void Slip::update_rho_hard(double *statev, double dtime, double temperature){
    // Calculate the coplanar interaction in FCC 111 slip systems
    double d_term_coplanar = 0, coeff_coplanar = harden_params[15];
    int index_1 = 1000, index_2 = 1000;
    for (int mode_id = 0; mode_id < total_mode_num; ++mode_id){
        if (mode_sys[mode_id]->type != slip) continue;
        if (mode_sys[mode_id]->num == num) continue;
        int inter_mode = interaction_mat(num, mode_id);
        if (inter_mode != 2) continue;
        if (index_1 == 1000) index_1 = mode_id;
        else index_2 = mode_id;
        if (index_1 != 1000 && index_2 != 1000) break;
    }
    double rho_1 = statev[sdv_ind(index_1, "DD")], rho_2 = statev[sdv_ind(index_2, "DD")], \
           rho_x = statev[sdv_ind(num, "DD")];
    double vel_1 = statev[sdv_ind(index_1, "SSR")] / (rho_1 * mode_sys[index_1]->burgers_vec.norm() * 1e-10), \
           vel_2 = statev[sdv_ind(index_2, "SSR")] / (rho_2 * mode_sys[index_2]->burgers_vec.norm() * 1e-10), \
           vel_x = statev[sdv_ind(num, "SSR")] / (rho_x * burgers_vec.norm() * 1e-10);
    double plus_term = rho_1 * vel_1 * sqrt(rho_2) + rho_2 * vel_2 * sqrt(rho_1);
    double minus_term = rho_1 * vel_1 * sqrt(rho_x) + rho_x * vel_x * sqrt(rho_1);
    minus_term += rho_2 * vel_2 * sqrt(rho_x) + rho_x * vel_x * sqrt(rho_2);
    d_term_coplanar = plus_term - minus_term;
    d_term_coplanar = d_term_coplanar * coeff_coplanar * dtime;
    if (d_term_coplanar + rho_x < 0) d_term_coplanar = -rho_x * 0.5;
    /* custom_var = d_term_coplanar; */
    statev[sdv_ind(num, "DD")] = rho_x + d_term_coplanar;
    rho_init = min(rho_x+d_term_coplanar, rho_init);
}
