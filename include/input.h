#include "sxcpp.h"
#include "pmode.h"
#include "slip.h"

inline Matrix6d read_elastic(const char* filename) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        throw std::runtime_error("Error opening file");
    }

    int row_num = 0, temp_idx = 0;
    double temp[6] = {0, 0, 0, 0, 0, 0};
    char temp_str[256];
    Matrix6d modulus = Matrix6d::Zero();

    while (row_num < 6 && fgets(temp_str, sizeof(temp_str), file)) {
        temp_idx = 0;
        char* token = strtok(temp_str, " ");
        while (token != nullptr && temp_idx < 6) {
            temp[temp_idx++] = atof(token);
            token = strtok(nullptr, " ");
        }
        modulus.row(row_num++) << temp[0], temp[1], temp[2], temp[3], temp[4], temp[5];
    }

    fclose(file);
    return modulus;
}

inline Matrix3d read_lattice(const char* filename) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        throw std::runtime_error("Error opening file");
    }

    char temp_str[256];
    // Skip the first 6 lines
    for (int i = 0; i < 6; ++i) {
        if (!fgets(temp_str, sizeof(temp_str), file)) {
            fclose(file);
            throw std::runtime_error("Error reading file: not enough lines");
        }
    }

    int row_num = 0, temp_idx = 0;
    double temp[3] = {0, 0, 0};
    Matrix3d lattice_mat = Matrix3d::Zero();

    // Read the next 3 lines to get the lattice matrix
    while (row_num < 3 && fgets(temp_str, sizeof(temp_str), file)) {
        temp_idx = 0;
        char* token = strtok(temp_str, " ");
        while (token != nullptr && temp_idx < 3) {
            temp[temp_idx++] = atof(token);
            token = strtok(nullptr, " ");
        }
        lattice_mat.row(row_num++) << temp[0], temp[1], temp[2];
    }

    fclose(file);
    return lattice_mat;
}

inline void read_pmodes(const char* filename) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        throw std::runtime_error("Error opening file");
    }

    char temp_str[256];
    double temp_info[MAX_MODE_NUM][6] = {0};
    double temp_harden[MAX_MODE_NUM] = {0};
    double temp_latent[MAX_LATENT_NUM] = {0};
    double temp_num = 0;
    int temp_idx = 0;
    mode_type temp_type;

    // Skip the first 9 lines
    for (int i = 0; i < 9; ++i) {
        if (!fgets(temp_str, sizeof(temp_str), file)) {
            fclose(file);
            throw std::runtime_error("Error reading file: not enough lines");
        }
    }

    while (total_mode_num != MAX_MODE_NUM) {
        // Read pmode type: slip, twin
        if (!fgets(temp_str, sizeof(temp_str), file)) break;
        if (strstr(temp_str, "Slip") || strstr(temp_str, "slip")) {
            temp_type = slip;
            flag_harden = 1;
        } else if (strstr(temp_str, "Voce") || strstr(temp_str, "voce")) {
            temp_type = slip;
            flag_harden = 0;
        } else if (strstr(temp_str, "Twin") || strstr(temp_str, "twin")) {
            temp_type = twin;
        } else {
            temp_type = undefined;
        }

        // Read numbers of pmode system
        if (!fgets(temp_str, sizeof(temp_str), file)) break;
        temp_num = atof(temp_str);

        // Read plane normal and burgers vector
        for (int mode_child_idx = 0; mode_child_idx < temp_num; ++mode_child_idx) {
            if (!fgets(temp_str, sizeof(temp_str), file)) break;
            temp_idx = 0;
            char* token = strtok(temp_str, " ");
            while (token != nullptr && temp_idx < 6) {
                temp_info[mode_child_idx][temp_idx++] = atof(token);
                token = strtok(nullptr, " ");
            }
        }

        // Read hardening parameters
        if (!fgets(temp_str, sizeof(temp_str), file)) break;
        temp_idx = 0;
        char* token = strtok(temp_str, " ");
        while (token != nullptr && temp_idx < MAX_HARDEN_NUM) {
            temp_harden[temp_idx++] = atof(token);
            token = strtok(nullptr, " ");
        }

        // Read latent parameters
        if (!fgets(temp_str, sizeof(temp_str), file)) break;
        temp_idx = 0;
        token = strtok(temp_str, " ");
        while (token != nullptr && temp_idx < MAX_LATENT_NUM) {
            temp_latent[temp_idx++] = atof(token);
            token = strtok(nullptr, " ");
        }

        // Create pmode system
        Vector6d slip_info = Vector6d::Zero();
        HardenVec hardens = HardenVec::Zero();
        LatentVec latents = LatentVec::Zero();
        for (int i = 0; i < MAX_HARDEN_NUM; ++i) hardens(i) = temp_harden[i];
        for (int i = 0; i < MAX_LATENT_NUM; ++i) latents(i) = temp_latent[i];

        for (int mode_child_idx = 0; mode_child_idx < temp_num; ++mode_child_idx) {
            int pmode_num = total_mode_num + mode_child_idx;
            for (int i = 0; i < 6; ++i) slip_info(i) = temp_info[mode_child_idx][i];
            if (temp_type == slip) {
                slip_array[pmode_num].initialize(pmode_num, slip_info, hardens, latents, lattice_vec, 1.0);
                mode_sys[pmode_num] = &slip_array[pmode_num];
            } else if (temp_type == twin) {
                slip_array[pmode_num].initialize(pmode_num, slip_info, hardens, latents, lattice_vec, 1.0);
                mode_sys[pmode_num] = &slip_array[pmode_num];
            } else {
                std::cerr << "Undefined mode type" << std::endl;
                fclose(file);
                throw std::runtime_error("SXCpp Error: Undefined mode type");
            }
        }
        total_mode_num += temp_num;
    }

    fclose(file);
}
