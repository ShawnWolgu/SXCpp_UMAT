#include "sxcpp.h"
#include "pmode.h"
#include "slip.h"

inline Matrix6d read_elastic(ifstream &is){
    int row_num = 0, temp_idx = 0;
    double temp[6] = {0,0,0,0,0,0};
    string temp_str;
    Matrix6d modulus;
    for(;row_num !=6; ++row_num){
        temp_idx = 0;
        getline(is, temp_str);
        stringstream stream(temp_str);
        while(!stream.eof() && temp_idx!=6) stream >> temp[temp_idx++];
        modulus.row(row_num) << temp[0], temp[1], temp[2], temp[3], temp[4], temp[5];
    }
    /* cout << "Read Elastic Modulus:" << endl << modulus << endl; */
    return modulus;
}

inline Matrix3d read_lattice(ifstream &is){
    int row_num = 0, temp_idx = 0;
    double temp[3] = {0,0,0};
    string temp_str;
    Matrix3d lattice_mat;
    for(;row_num !=3; ++row_num){
        temp_idx = 0;
        getline(is, temp_str);
        stringstream stream(temp_str);
        while(!stream.eof() && temp_idx!=3) stream >> temp[temp_idx++];
        lattice_mat.row(row_num) << temp[0], temp[1], temp[2];
    }
    /* cout << "Read Elastic Modulus:" << endl << modulus << endl; */
    return lattice_mat;
}

inline void read_pmodes(ifstream &is){
    string temp_str;
    double temp_num = 0;
    double temp_info[MAX_MODE_NUM][6] = {0};
    double temp_harden[MAX_MODE_NUM] = {0};
    double temp_latent[MAX_LATENT_NUM] = {0};
    mode_type temp_type;
    int temp_idx = 0;
    for(;total_mode_num!=MAX_MODE_NUM;){
        //Read pmode type: slip, twin
        getline(is, temp_str);
        if (is.eof()) break;
        if((temp_str.find("Slip") != temp_str.npos) || (temp_str.find("slip") != temp_str.npos)) {
            temp_type = slip;
            flag_harden = 1;
        } 
        else if((temp_str.find("Voce") != temp_str.npos) || (temp_str.find("voce") != temp_str.npos)) {
            temp_type = slip;
            flag_harden = 0;
        } 
        else if((temp_str.find("Twin") != temp_str.npos) || (temp_str.find("twin") != temp_str.npos)) {
            temp_type = twin;
        } 
        else {
            temp_type = undefined;
        }
        //Read numbers of pmode system
        getline(is, temp_str);
        stringstream stream(temp_str);
        stream >> temp_num;
        stream.clear();
        //Read plane normal and burgers vector
        for (int mode_child_idx = 0; mode_child_idx < temp_num; ++mode_child_idx){
            temp_idx = 0;
            getline(is, temp_str);
            stream.str(temp_str);
            while(!stream.eof() && temp_idx!=6) stream >> temp_info[mode_child_idx][temp_idx++];
            stream.clear();
            stream.str("");
        }
        //Read hardening parameters
        getline(is, temp_str);
        stream.str(temp_str);
        temp_idx = 0;
        while(!stream.eof() && temp_idx!=MAX_HARDEN_NUM) stream >> temp_harden[temp_idx++];
        stream.clear();
        //Read latent parameters
        getline(is, temp_str);
        stream.str(temp_str);
        temp_idx = 0;
        while(!stream.eof() && temp_idx!=MAX_LATENT_NUM) stream >> temp_latent[temp_idx++];
        //Create pmode system
        Vector6d slip_info = Vector6d::Zero();
        HardenVec hardens = HardenVec::Zero();
        LatentVec latents = LatentVec::Zero();
        for(int i=0; i<MAX_HARDEN_NUM; ++i) hardens(i) = temp_harden[i];
        for(int i=0; i<MAX_LATENT_NUM; ++i) latents(i) = temp_latent[i];
        for (int mode_child_idx = 0; mode_child_idx < temp_num; ++mode_child_idx){
            int pmode_num = total_mode_num + mode_child_idx;
            for(int i=0; i<6; ++i) slip_info(i) = temp_info[mode_child_idx][i];
            if (temp_type == slip) 
                mode_sys[pmode_num] = make_unique<Slip>(pmode_num, slip_info, hardens, latents, lattice_vec, 1.0);
            else if (temp_type == twin) 
                mode_sys[pmode_num] = make_unique<Slip>(pmode_num, slip_info, hardens, latents, lattice_vec, 1.0);
                /* mode_sys[pmode_num] = Twin(pmode_num, slip_info, hardens, latents, lattice_vec, f_active); */
            else {
                cout << "Undefined mode type" << endl;
                break;
            }
        }
        total_mode_num += temp_num;
    }
}
