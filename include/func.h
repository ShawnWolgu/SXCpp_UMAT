#pragma once

#include "sxcpp.h"

Matrix6d cal_rotation_trans_6d_for_stiff(Matrix3d M);
Matrix6d cal_rotation_trans_6d_for_compl(Matrix3d M);
int get_interaction_mode(Vector3d burgers_i, Vector3d plane_i, Vector3d burgers_j, Vector3d plane_j);

inline Matrix6d Jacobian_Matrix(double* ddsdde){
    Matrix6d C;
    for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < 6; j++)
        {
            C(i, j) = ddsdde[i * 6 + j];
        }
    }
    return C;
}

inline void ddsdde_from_matrix(Matrix<double, 6, 6> Jacobian, double* ddsdde){
    for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < 6; j++)
        {
            ddsdde[i * 6 + j] = Jacobian(i, j);
        }
    }
}

inline Vector6d change_basis_order(Vector6d C){
    Vector6d C_out;
    Matrix6d change_matrix;
    change_matrix << 1, 0, 0, 0, 0, 0,
                     0, 1, 0, 0, 0, 0,
                     0, 0, 1, 0, 0, 0,
                     0, 0, 0, 0, 0, 1,
                     0, 0, 0, 0, 1, 0,
                     0, 0, 0, 1, 0, 0;
    C_out = change_matrix * C;
    return C_out;
}

inline Matrix6d change_basis_order(Matrix6d C){
    Matrix6d C_out;
    Matrix6d change_matrix;
    change_matrix << 1, 0, 0, 0, 0, 0,
                     0, 1, 0, 0, 0, 0,
                     0, 0, 1, 0, 0, 0,
                     0, 0, 0, 0, 0, 1,
                     0, 0, 0, 0, 1, 0,
                     0, 0, 0, 1, 0, 0;
    C_out = change_matrix * C * change_matrix.transpose();
    return C_out;
}

inline Matrix3d Euler_trans(double phi1, double Phi, double phi2){
    Matrix3d Mout;
    phi1 = phi1 * pi / 180; Phi = Phi * pi / 180; phi2 = phi2 * pi / 180;
    double SPH,STH,STM;
    double CPH,CTH,CTM;
    SPH=sin(phi1); CPH=cos(phi1);
    STH=sin(Phi); CTH=cos(Phi);
    STM=sin(phi2); CTM=cos(phi2);
    Mout(0,0)=CTM*CPH-SPH*STM*CTH;
    Mout(1,0)=-STM*CPH-SPH*CTM*CTH;
    Mout(2,0)=SPH*STH;
    Mout(0,1)=CTM*SPH+CPH*STM*CTH;
    Mout(1,1)=-SPH*STM+CPH*CTM*CTH;
    Mout(2,1)=-STH*CPH;
    Mout(0,2)=STH*STM;
    Mout(1,2)=CTM*STH;
    Mout(2,2)=CTH;
    return Mout;
}

inline Vector3d Euler_trans(Matrix3d euler_matrix){
    Vector3d Vout;
    double PH,TH,TM;
    double STH;

    TH = acos(euler_matrix(2,2));
    if(abs(euler_matrix(2,2)) >= 0.99)
    {
        TM = 0;
        PH = atan2(euler_matrix(0,1),euler_matrix(0,0));
    }
    else
    {
        STH = sin(TH);
        TM = atan2(euler_matrix(0,2)/STH,euler_matrix(1,2)/STH);
        PH = atan2(euler_matrix(2,0)/STH,-euler_matrix(2,1)/STH); 
    }
    Vout << PH,TH,TM;
    Vout = Vout * 180 / pi; // in degree
    return Vout;
}

inline Matrix6d rotate_6d_stiff_modu(Matrix6d modulus, Matrix3d rotate_matrix){
    Matrix6d M66 = cal_rotation_trans_6d_for_stiff(rotate_matrix);
    return M66*modulus*M66.transpose();
}

inline Matrix6d cal_rotation_trans_6d_for_stiff(Matrix3d M){
   Matrix6d M66;
    double xx,yy,zz,yz,xz,xy,yx,zx,zy;
    xx = M(0,0); 
    yy = M(1,1);
    zz = M(2,2);
    yz = M(1,2);
    xz = M(0,2);
    xy = M(0,1);
	yx = M(1,0);
	zx = M(2,0);
	zy = M(2,1);

    for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
    M66(i,j) = M(i,j)*M(i,j);

    M66(0,3) = 2*xy*xz;
    M66(0,4) = 2*xx*xz;
    M66(0,5) = 2*xx*xy;

    M66(1,3) = 2*yy*yz;
    M66(1,4) = 2*yx*yz;
    M66(1,5) = 2*yx*yy;

    M66(2,3) = 2*zy*zz;
    M66(2,4) = 2*zx*zz;
    M66(2,5) = 2*zx*zy;

    M66(3,0) = yx*zx;
    M66(3,1) = yy*zy;
    M66(3,2) = yz*zz;

    M66(4,0) = xx*zx;
    M66(4,1) = xy*zy;
    M66(4,2) = xz*zz;

    M66(5,0) = xx*yx;
    M66(5,1) = xy*yy;
    M66(5,2) = xz*yz;

    M66(3,3) = yy*zz+yz*zy;
    M66(3,4) = yx*zz+yz*zx;
    M66(3,5) = yx*zy+yy*zx;

    M66(4,3) = xy*zz+xz*zy;
    M66(4,4) = xx*zz+xz*zx;
    M66(4,5) = xx*zy+xy*zx;

    M66(5,3) = xy*yz+xz*yy;
    M66(5,4) = xx*yz+xz*yx;
    M66(5,5) = xx*yy+xy*yx;

	return M66;
}

inline Matrix6d rotate_6d_compl_modu(Matrix6d modulus, Matrix3d rotate_matrix)
{
    Matrix6d M66 = cal_rotation_trans_6d_for_compl(rotate_matrix);
    return M66*modulus*M66.transpose();
}

inline Matrix6d cal_rotation_trans_6d_for_compl(Matrix3d M){
	Matrix6d M66;
    double xx,yy,zz,yz,xz,xy,yx,zx,zy;
    xx = M(0,0); 
    yy = M(1,1);
    zz = M(2,2);
    yz = M(1,2);
    xz = M(0,2);
    xy = M(0,1);
	yx = M(1,0);
	zx = M(2,0);
	zy = M(2,1);

    for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
            M66(i,j) = M(i,j)*M(i,j);

    M66(0,3) = xy*xz;
    M66(0,4) = xx*xz;
    M66(0,5) = xx*xy;

    M66(1,3) = yy*yz;
    M66(1,4) = yx*yz;
    M66(1,5) = yx*yy;

    M66(2,3) = zy*zz;
    M66(2,4) = zx*zz;
    M66(2,5) = zx*zy;

    M66(3,0) = 2*yx*zx;
    M66(3,1) = 2*yy*zy;
    M66(3,2) = 2*yz*zz;

    M66(4,0) = 2*xx*zx;
    M66(4,1) = 2*xy*zy;
    M66(4,2) = 2*xz*zz;

    M66(5,0) = 2*xx*yx;
    M66(5,1) = 2*xy*yy;
    M66(5,2) = 2*xz*yz;

    M66(3,3) = yy*zz+yz*zy;
    M66(3,4) = yx*zz+yz*zx;
    M66(3,5) = yx*zy+yy*zx;

    M66(4,3) = xy*zz+xz*zy;
    M66(4,4) = xx*zz+xz*zx;
    M66(4,5) = xx*zy+xy*zx;

    M66(5,3) = xy*yz+xz*yy;
    M66(5,4) = xx*yz+xz*yx;
    M66(5,5) = xx*yy+xy*yx;

	return M66;
}

inline void processPath(char* path, const char* terminationSequence) {
    const size_t terminationLength = strlen(terminationSequence);
    const size_t maxPathLength = 200;
    // Find the length of the valid path part (non-empty part)
    size_t len = 0;
    for (size_t i = 0; i < maxPathLength && path[i] != '\0'; ++i) {
        if (path[i] != ' ') {
            len = i + 1;
        }
    }

    // Add the termination sequence after the valid path
    if (len + terminationLength < 200) { // Ensure there's enough space for the sequence and null terminator
        strcpy(path + len, terminationSequence);
    } else {
        std::cerr << "Error: Not enough space to add termination sequence." << std::endl;
    }
}

inline Vector6d tensor_trans_order(Matrix3d tensor){
    /** 
     * 3*3 Matrix --> 1*6 Vector: [00,11,22,12,02,01] 
     */
    return Vector6d{{tensor(0,0), tensor(1,1), tensor(2,2), tensor(1,2), tensor(0,2), tensor(0,1)}};
}

inline Matrix<double,9,1> tensor_trans_order_9(Matrix3d tensor){
    /** 
     * 3*3 Matrix --> 9*1 Vector: [00,11,22,12,02,01,21,20,10] 
     */
    return Matrix<double,9,1> {{tensor(0,0), tensor(1,1), tensor(2,2), tensor(1,2), tensor(0,2), tensor(0,1), tensor(2,1), tensor(2,0), tensor(1,0)}};
}

inline Matrix3d tensor_trans_order(Vector6d tensor){
    /** 
     * 1*6 Vector --> 3*3 Matrix: [[0,5,4],[5,1,3],[4,3,2]] 
     */
    return Matrix3d{{tensor(0), tensor(5), tensor(4)}, {tensor(5), tensor(1), tensor(3)}, {tensor(4), tensor(3), tensor(2)}};
}

inline Matrix3d tensor_trans_order_9(Matrix<double,9,1> tensor){
    /** 
     * 9*1 Vector --> 3*3 Matrix: [[0,5,4],[8,1,3],[7,6,2]] 
     */
    return Matrix3d{{tensor(0), tensor(5), tensor(4)}, {tensor(8), tensor(1), tensor(3)}, {tensor(7), tensor(6), tensor(2)}};
}

inline Matrix6d get_C_ij_pri(Matrix6d elastic_modulus, Vector6d &stress_6d){
    Matrix6d C_ij_pri;
    C_ij_pri << stress_6d, stress_6d, stress_6d, Vector6d::Zero(), Vector6d::Zero(), Vector6d::Zero();
    C_ij_pri = elastic_modulus - C_ij_pri;
    return C_ij_pri;
}

inline Matrix3d Rodrigues(Matrix3d spin_elas){
    Matrix3d exp_skew;
    Vector3d omegas;
    omegas << spin_elas(2,1), spin_elas(0,2), spin_elas(1,0);
    double omega_norm = omegas.norm();
    if(omega_norm < 1e-30)
        exp_skew = Matrix3d::Identity();
    else
    {
        Matrix3d term_1 = spin_elas/omega_norm * sin(omega_norm);
        Matrix3d term_2 = (spin_elas/omega_norm) * (spin_elas/omega_norm) * (1.0 - cos(omega_norm));
        exp_skew =  Matrix3d::Identity()+ term_1 + term_2;
    }
    return exp_skew;
}

inline Matrix3d tensor_rot_to_CryCoord(Matrix3d tensor, Matrix3d orientation){
    return orientation * tensor * orientation.transpose();
}

inline Matrix3d tensor_rot_to_RefCoord(Matrix3d tensor, Matrix3d orientation){
    return orientation.transpose() * tensor * orientation;
}

inline int sign(double x){
    return (0 < x) - (x < 0);
}

inline int heaviside(double x){
    return (x > 0);
}

inline double set_precision(double num, int prec){
    if (num == 0) return 0.0;
    double expo = floor(log10(abs(num)));
    double result = sign(num) * round(abs(num) * pow(10,-expo + prec -1))/pow(10,-expo + prec -1);
    return result;
}

inline double calc_equivalent_value(Matrix3d mat){
    Matrix3d dev_mat = mat - Matrix3d::Identity() * mat.trace();
    double sum = (dev_mat.cwiseProduct(dev_mat)).sum();
    return sqrt(2./3. * sum);
}

inline double calc_first_principal(Matrix3d mat){
    Vector3d eigens = mat.eigenvalues().cwiseAbs();
    return eigens.maxCoeff();
}

inline Vector6d set_precision(Vector6d &num, int prec){
    Vector6d result; result << num;
    for (auto &inum : result){
        if (inum != 0){
            double expo = floor(log10(abs(inum)));
            inum = sign(inum) * round(abs(inum) * pow(10,-expo + prec -1))/pow(10,-expo + prec -1);
        }
    }
    return result;
}

inline void cut_precision(Matrix3d &mat, int prec){
    double ref_value = mat.cwiseAbs().maxCoeff();
    if (ref_value != 0){
	double expo = floor(log10(ref_value)) - prec;
	mat = round(mat.array() * pow(10,-expo)).matrix() * pow(10,expo);
    }
}

inline double calc_relative_error(Vector6d &v1, Vector6d &v2){
    Vector6d v_error = (v1 - v2).cwiseAbs();
    Vector6d inv_v1 = (0.5 * v1 + 0.5 * v2);
    for (auto &iv : inv_v1) if (iv != 0) iv = 1/iv;
    double result = (v_error.cwiseProduct(inv_v1)).norm();
    return result;
}

inline double calc_relative_error(double x, double y){
    if (x == 0 && y == 0) return 0;
    else return abs(x)-abs(y)/abs(x);
}

inline double relative_std(VectorXd &vec){
    double mean = vec.mean();
    double sum = 0;
    for (int i = 0; i < vec.size(); i++) sum += pow(vec(i) - mean,2);
    return sqrt(sum/vec.size())/mean;
}

inline Vector3d get_plane_norm(Vector3d &plane_norm_disp, Matrix3d &lattice_vec){
    double h = plane_norm_disp(0), k = plane_norm_disp(1), l = plane_norm_disp(2);
    if (h == 0 && k == 0 && l == 0) return Vector3d::Zero();
    Vector3d h_vec = (h==0)? Vector3d(1e10,0,0) : Vector3d(1/h,0,0);
    Vector3d k_vec = (k==0)? Vector3d(0,1e10,0) : Vector3d(0,1/k,0);
    Vector3d l_vec = (l==0)? Vector3d(0,0,1e10) : Vector3d(0,0,1/l);
    Vector3d vec_one = ((k_vec - h_vec).transpose() * lattice_vec);
    Vector3d vec_two = ((l_vec - k_vec).transpose() * lattice_vec);
    Vector3d norm = vec_one.cross(vec_two);
    return norm / norm.norm();
}

inline double cal_cosine(Vector3d vec_i, Vector3d vec_j){
    return vec_i.dot(vec_j)/(vec_i.norm() * vec_j.norm());
}

/*
 * Return the dislocation interaction mode code between two slip system.
 * 0: No Junction, 1: Hirth Lock, 2: Coplanar Junction, 3: Glissile Junction, 4: Sessile Junction
 */
inline int get_interaction_mode(Vector3d burgers_i, Vector3d plane_i, Vector3d burgers_j, Vector3d plane_j){
    double perp = 0.02, prll = 0.98;
    double cos_b_angle = cal_cosine(burgers_i, burgers_j);
    if(abs(cos_b_angle) < perp) return 1;
    else {
        if(abs(cos_b_angle) > prll) return 0;
        else{
            if (abs(cal_cosine(plane_i, plane_j)) > prll) return 2;
            else{
                bool if_glide_i = (abs(cal_cosine(plane_i, burgers_i+burgers_j)) < perp);
                bool if_glide_j = (abs(cal_cosine(plane_j, burgers_i+burgers_j)) < perp);
                if (if_glide_i || if_glide_j) return 3;
                else return 4;
            }
        }
    }
}

inline int sdv_ind(int slip_num, string sdv_name){
    int start = 14 + slip_num * 5;
    //start+0: SSR; start+1: slope; start+2:ACC/DD/TVF;
    //start+3: CRSS/tauf; start+4: custom;
    if (sdv_name == "SSR") return start;
    else if (sdv_name == "slope") return start + 1;
    else if (sdv_name == "ACC") return start + 2;
    else if (sdv_name == "DD") return start + 2;
    else if (sdv_name == "TVF") return start + 2;
    else if (sdv_name == "CRSS") return start + 3;
    else if (sdv_name == "tauf") return start + 3;
    else if (sdv_name == "custom") return start + 4;
    else if (sdv_name == "mfp") return start + 4;
    else return -1;
}
