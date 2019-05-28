// Update 181225
#include "hamiltonian.hpp"
#include "config.hpp"
#include<stdlib.h>
#include<iostream>
#include<float.h>

using namespace std;

double dig2rad(double dig){ return M_PI*dig/180.0; }
double rad2dig(double rad){ return rad*180.0/M_PI; }
complex<double> SphericalHarmonicFunction(double x,double y,double z,int m){
    double r = sqrt(x*x + y*y + z*z);
    double theta = acos(z/r), phi = atan2(y,x);
    const double cY_0 = sqrt(5.0/(16.0*M_PI));
    const double cY_1 = sqrt(15.0/(8.0*M_PI));
    const double cY_2 = sqrt(15.0/(32.0*M_PI));
    
    if (r < DBL_EPSILON)
        return 0.0;
    else {
        switch(m){
            case -2: return cY_2*sin(theta)*sin(theta)*exp(-2*phi*I);
            case -1: return cY_1*(+1)*sin(theta)*cos(theta)*exp(-phi*I);
            case  0: return cY_0*(3*cos(theta)*cos(theta)-1);
            case  1: return cY_1*(-1)*sin(theta)*cos(theta)*exp(phi*I);
            case  2: return cY_2*sin(theta)*sin(theta)*exp(2*phi*I);
            default:
                cout << "Wrong index of Spherical Function.";
                return 0;
            
        }
    }
}


// MonoLayer
MonoLayer::MonoLayer(){
    t_ip=1.0;
    t_op=1.0;
    ex_sd=1.0;
    soc_r=1.0;
    rad_p=0.0;
    rad_a=0.0;
}
const MonoLayer MonoLayer::operator *(const MonoLayer& rhs) const{
    MonoLayer tmp;
    tmp.t_ip = t_ip * rhs.t_ip;
    tmp.t_op = t_op * rhs.t_op;
    tmp.ex_sd = ex_sd * rhs.ex_sd;
    tmp.soc_r = soc_r * rhs.soc_r;
    tmp.rad_p = rad_p + rhs.rad_p;
    tmp.rad_a = rad_a + rhs.rad_a;
    return tmp;
}

MonoLayer& MonoLayer::operator *=(const MonoLayer& rhs){
    MonoLayer tmp;
    t_ip *= rhs.t_ip;
    t_op *= rhs.t_op;
    ex_sd *= rhs.ex_sd;
    soc_r *= rhs.soc_r;
    rad_p += rhs.rad_p;
    rad_a += rhs.rad_a;
    return *this;
}

// D_orbital Impurity
Impurity::Impurity(){
    d_orbital_level = 0;
    d_orbital_width = 0;
    d_orbital_exchange = 0;
    d_orbital_soi = 0;
    d_orbital_screened_charge = 0;
    d_orbital_tetra = 0;
    d_orbital_cf = 0;
}

SpinRotation::SpinRotation(){ SpinRotation(0.0,0.0); }
SpinRotation::SpinRotation(double theta, double phi){ SpinRotation('z',theta,phi); }
SpinRotation::SpinRotation(char axis, double theta, double phi){
    if(axis == 'z'){
        rotmat[0][0] = cos(theta);
        rotmat[0][1] = sin(theta)*exp(-I*phi);
        rotmat[1][0] = sin(theta)*exp(I*phi);
        rotmat[1][1] = -cos(theta);
    }
    else if (axis == 'x')
        for (int i=0;i<2;i++)
            for(int j=0;j<2;j++)
                rotmat[i][j] = cos(theta)*Pauli[1][i][j] 
                                + sin(theta)*cos(phi)*Pauli[2][i][j] 
                                + sin(theta)*sin(phi)*Pauli[3][i][j];
    else
        for (int i=0;i<2;i++)
            for(int j=0;j<2;j++)
                rotmat[i][j] = sin(theta)*sin(phi)*Pauli[1][i][j] 
                                + cos(theta)*Pauli[2][i][j] 
                                + sin(theta)*cos(phi)*Pauli[3][i][j];
}
complex<double> SpinRotation::getRot(int i,int j){ return rotmat[i][j]; }


// Deffinition of Hamiltonian
Hamiltonian::Hamiltonian(){
    Initialize();
}

void Hamiltonian::Initialize() {
    root.t_ip = INIT_IPHOPPING;
    root.ex_sd = INIT_EXCHANGE;
    root.soc_r = INIT_SOC;
    root.rad_p = dig2rad(INIT_THETA);
    root.rad_a = dig2rad(INIT_PHI);
    root.t_op = INIT_OPHOPPING;
    
    num_layer = INIT_LAYER;
    gamma = INIT_GAMMA;
    epsf = INIT_FERMI_ENERGY;
    filling = INIT_FILLING;
    
    
    imp.d_orbital_level = INIT_IMP_LEVEL;
    imp.d_orbital_width = INIT_IMP_WIDTH;
    imp.d_orbital_exchange = INIT_IMP_EXCHANGE;
    imp.d_orbital_soi = INIT_IMP_SOI;
    //~ imp.d_orbital_screened_charge = INIT_IMP_CHARGE;
    imp.d_orbital_cf = INIT_IMP_CF;
    imp.d_orbital_tetra = INIT_IMP_TETRA_CF;
    
    n_imp = INIT_NIMP;
    v_sd = INIT_HYB_SD;
    
    Hd.SetMatrix(10);
    Gd_r.SetMatrix(10);
    Gd_a.SetMatrix(10);
}

// -- 00:Exchange,          01:Rashba-SOI,      02:In-Plane HOPPING,    03:Out-of-Plene Hopping,
// -- 04:Polar-angle(theta):[0:PI],       05:Azmual-angle(phi):[0:2PI], 06:4s Lorentzian width
// -- 10:Molecularfield,    11:Atomic-SOI,      12:d-orbital level,     13:Cubic cristal field splitting 
// -- 14:3d Lorentzian width 15:Tetragonal crystal field splitting
// -- 20:Impurity density,  21:Hybridization, 
// -- 90:Fermi energy(Chemical potential),     91:Filling(particle/unit cell):[0.0:2.0],    92:Layer [1:99]
char Hamiltonian::setParam(double a){
    char update_flag=1;
    switch (INPUT) {
            
        case 00:
            root.ex_sd = a;
            break;
        case 01:
            root.soc_r = a;
            break;
        case 02:
            root.t_ip = a;
            break;
        case 03:
            root.t_op = a;
            break;
        case 04:
            root.rad_p = dig2rad(a);
            break;
        case 05:
            root.rad_a = dig2rad(a);
            break;

        case 10:
            imp.d_orbital_exchange = a;
            break;
        case 11:
            imp.d_orbital_soi = a;
            break;
        case 12:
            imp.d_orbital_level = a;
            break;
        case 13:
            imp.d_orbital_cf = a;
            break;
        case 14:
            imp.d_orbital_width = a;
            break;
        case 15:
            imp.d_orbital_tetra = a;
            break;
        
        case 20:
            n_imp = a;
            break;
        case 21:
            v_sd = a;
            break;

        case 90:
            epsf = a;
            update_flag=0;
            break;
        case 91:
            filling = a;
            update_flag=0;
            break;
        case 92:
            num_layer = (int)a;
            break;
    }
    
    setSD_Impurity();
    return update_flag;
}

    
void Hamiltonian::printParam(){
    printParam('I');
    printParam('E');
}

void Hamiltonian::printParam(char option){
    if (option == 'I'){
        cout << "\nIntrinsic Parameters..." <<endl;
        cout << "In-plane hopping: " << root.t_ip << endl;
        cout << "Out-of-plane hopping: " << root.t_op << endl;
        cout << "Exchange splitting: " << root.ex_sd << endl;
        cout << "Rashba-SOI coupling constant: " << root.soc_r << endl;
        cout << "Magnetizatoin orientetion: (" << root.rad_p << "," << root.rad_a << ")" << endl;
    }
    else if (option == 'E'){
        cout << "\nExtrinsic Parameters..." <<endl;
        cout << "Number of Layer: " << num_layer << endl;
        cout << "Spectrum width: " << gamma << endl;
        cout << "Band filling: " << filling << endl;
        cout << "Fermi-energy: " << epsf << endl;
    }
}

void Hamiltonian::setDirection(char dir){
    if (dir == 'x') setDirection(90.0,0.0);
    else if (dir == 'y')setDirection(90.0,90.0);
    else if (dir == 'z')setDirection(0.0,0.0);
}
void Hamiltonian::setDirection(double theta, double phi){
    root.rad_p = dig2rad(theta);
    if(NODAL_CONDUCTION==1)root.rad_a += dig2rad(phi + 45.0);
    else root.rad_a = dig2rad(phi);
    setSD_Impurity();
}

void Hamiltonian::setFermiEnergy(double x){ epsf = x; }
int Hamiltonian::getLayer(){ return num_layer; }
MonoLayer Hamiltonian::getRootParam(){ return root; }
double Hamiltonian::getLifetime(){ return gamma; }
double Hamiltonian::getFermiEnergy(){ return epsf; }
double Hamiltonian::getFilling(){ return filling; }
char Hamiltonian::check_Integral_type(){ return integral_type; }
string Hamiltonian::H_Name(){ return h_name; }

void Hamiltonian::OutputHamiltonian(ofstream& fout){
    OutputStructure(fout); 
    
    fout << "# Impurity Parameter" << endl;
    fout << "# n_imp, Ed(t2g), EX_Hund, ASOI, Cubic_CF, Tetra_CF, gamma, |Vsd|" << endl;
    fout << "# " << n_imp << ", ";
    fout << imp.d_orbital_level << ", ";
    fout << imp.d_orbital_exchange << ", ";
    fout << imp.d_orbital_soi << ", ";
    fout << imp.d_orbital_cf << ", ";
    fout << imp.d_orbital_tetra << ",";
    fout << imp.d_orbital_width << ",";
    fout << v_sd << endl;
    fout << endl;
    
    fout << "# Environment Parameter" << endl;
    fout << "# EF, filling" << endl;
    fout << "# " << INIT_FERMI_ENERGY << ", ";
    fout << INIT_FILLING << endl;
    fout << endl;
    
    fout << "# Numeriacl Parameter" << endl;
    fout << "# k_min/PI, k_max/PI, k_div, E_min, E_max, E_div, const_EF, Omit_sea_term"  << endl;
    fout << "# " << K_MIN/M_PI << ", ";
    fout << K_MAX/M_PI << ", ";
    fout << K_DIV << ", ";
    fout << E_MIN << ", ";
    fout << E_MAX << ", ";
    fout << E_DIV << ", ";
    fout << EF_HALFWAY_MODE << ", ";
    fout << SF_HALFWAY_MODE << endl;
    fout << endl;
}



// Impurity
void Hamiltonian::setSD_Impurity(){
    double theta = root.rad_p, phi = root.rad_a;
    double xi_2 = 0.5*imp.d_orbital_soi, delta = imp.d_orbital_exchange, Et = imp.d_orbital_level, DE = 0.5*imp.d_orbital_cf,Dtetra=imp.d_orbital_tetra;
    SpinRotation R(ROTATION_AXIS,theta,phi);
    
    // Modulate imp. level as bottom level of d_epsilon.
    if( imp.d_orbital_cf < 0) Et -= imp.d_orbital_cf;
    
    Hd.SetMatrix(10);
#ifndef PERTUB
    for(int i=0;i<5;i++){
        for(int j=0;j<5;j++){
            int m = i-2, l = j-2;
            
            if(m == l){
                Hd.Element(2*i,2*j) = Et -delta*R.getRot(0,0) +xi_2 * l;
                Hd.Element(2*i+1,2*j+1) = Et -delta*R.getRot(1,1) -xi_2 * l;
                Hd.Element(2*i,2*j+1) = -delta*R.getRot(0,1);
                Hd.Element(2*i+1,2*j) = -delta*R.getRot(1,0);
                if ( m == 0 ){
                    Hd.Element(2*i,2*j) += (2*DE+Dtetra);
                    Hd.Element(2*i+1,2*j+1) += (2*DE+Dtetra);
                }
                if ( (m == 1) || (m == -1) ){
                    Hd.Element(2*i,2*j) += Dtetra;
                    Hd.Element(2*i+1,2*j+1) += Dtetra;
                }
            }
            else if (m == l+1){
                Hd.Element(2*i,2*j) = 0.0;
                Hd.Element(2*i+1,2*j+1) = 0.0;
                Hd.Element(2*i+1,2*j) = xi_2 * sqrt((2-l)*(2+l+1));
                Hd.Element(2*i,2*j+1) = 0.0;
            }
            else if (m == l-1){
                Hd.Element(2*i,2*j) = 0.0;
                Hd.Element(2*i+1,2*j+1) = 0.0;
                Hd.Element(2*i+1,2*j) = 0.0;
                Hd.Element(2*i,2*j+1) = xi_2 * sqrt((2+l)*(2-l+1));
            }
            else {
                Hd.Element(2*i,2*j) = 0.0;
                Hd.Element(2*i+1,2*j+1) = 0.0;
                Hd.Element(2*i,2*j+1) = 0.0;
                Hd.Element(2*i+1,2*j) = 0.0;
            }
            
            if ( (m*l == 4) || (m*l == -4) ){
                Hd.Element(2*i,2*j) += DE;
                Hd.Element(2*i+1,2*j+1) += DE;
            }
        }
    }
    
    //Hd.PrintMatrix();
    
    
#else    
     // Non-perturbation Hamiltonian
     double A=DE*0;
    for(int i=0;i<5;i++){
        for(int j=0;j<5;j++){
            int m = i-2, l = j-2;
            if(m == l){
                Hd.Element(2*i,2*j) = Et -delta*R.getRot(0,0);
                Hd.Element(2*i+1,2*j+1) = Et -delta*R.getRot(1,1);
                Hd.Element(2*i,2*j+1) = -delta*R.getRot(0,1);
                Hd.Element(2*i+1,2*j) = -delta*R.getRot(1,0);
                if ( m == 0 ){
                    Hd.Element(2*i,2*j) += 2*A;
                    Hd.Element(2*i+1,2*j+1) += 2*A;
                }
            }
            else if (m == l+1){
                Hd.Element(2*i,2*j) = 0.0;
                Hd.Element(2*i+1,2*j+1) = 0.0;
                Hd.Element(2*i+1,2*j) = 0.0;
                Hd.Element(2*i,2*j+1) = 0.0;
            }
            else if (m == l-1){
                Hd.Element(2*i,2*j) = 0.0;
                Hd.Element(2*i+1,2*j+1) = 0.0;
                Hd.Element(2*i+1,2*j) = 0.0;
                Hd.Element(2*i,2*j+1) = 0.0;
            }
            else {
                Hd.Element(2*i,2*j) = 0.0;
                Hd.Element(2*i+1,2*j+1) = 0.0;
                Hd.Element(2*i,2*j+1) = 0.0;
                Hd.Element(2*i+1,2*j) = 0.0;
            }
            if ( (m*l == 4) || (m*l == -4) ){
                Hd.Element(2*i,2*j) += A;
                Hd.Element(2*i+1,2*j+1) += A;
            }
        }
    }
#endif    
}

void Hamiltonian::setD_GreenFunc(double eps){
    double wid = imp.d_orbital_width;
    
    Gd_r.SetMatrix(10);
    Gd_a.SetMatrix(10);
#ifndef PERTUB
    for(int i=0;i<10;i++)
        for(int j=0;j<10;j++){
            if (i==j){
                Gd_r.Element(i,j) = eps - Hd.Element(i,j) + I*wid;
                Gd_a.Element(i,j) = eps - Hd.Element(i,j) - I*wid;
            }
            else {
                Gd_r.Element(i,j) = -Hd.Element(i,j);
                Gd_a.Element(i,j) = -Hd.Element(i,j);
            }
        }
    Gd_r.Inverse();
    Gd_a.Inverse();
    
    
    
#else    
#warning    // Peturbation Hamiltonian
    Complex_SquareMatrix K(10),L(10),M(10),N(10),O(10),P(10);
    for(int i=0;i<10;i++)
        for(int j=0;j<10;j++){
            if (i==j){
                K.Element(i,j) = eps - Hd.Element(i,j) + I*wid;
                L.Element(i,j) = eps - Hd.Element(i,j) - I*wid;
            }
            else {
                K.Element(i,j) = -Hd.Element(i,j);
                L.Element(i,j) = -Hd.Element(i,j);
            }
        }
    K.Inverse();
    L.Inverse();
    Complex_SquareMatrix HP(10);
    double xi_2 = 0.5*imp.d_orbital_soi,DE = 0.5*imp.d_orbital_cf;
    for(int i=0;i<5;i++){
        for(int j=0;j<5;j++){
            int m = i-2, l = j-2;
            
            if(m == l){
                HP.Element(2*i,2*j) = xi_2 * l;
                HP.Element(2*i+1,2*j+1) = -xi_2 * l;
                HP.Element(2*i,2*j+1) = 0.0;
                HP.Element(2*i+1,2*j) = 0.0;
                if ( m == 0 ){
                    HP.Element(2*i,2*j) += 2*DE;
                    HP.Element(2*i+1,2*j+1) += 2*DE;
                }
            }
            else if (m == l+1){
                HP.Element(2*i,2*j) = 0.0;
                HP.Element(2*i+1,2*j+1) = 0.0;
                HP.Element(2*i+1,2*j) = xi_2 * sqrt((2-l)*(2+l+1));
                HP.Element(2*i,2*j+1) = 0.0;
            }
            else if (m == l-1){
                HP.Element(2*i,2*j) = 0.0;
                HP.Element(2*i+1,2*j+1) = 0.0;
                HP.Element(2*i+1,2*j) = 0.0;
                HP.Element(2*i,2*j+1) = xi_2 * sqrt((2+l)*(2-l+1));
            }
            else {
                HP.Element(2*i,2*j) = 0.0;
                HP.Element(2*i+1,2*j+1) = 0.0;
                HP.Element(2*i,2*j+1) = 0.0;
                HP.Element(2*i+1,2*j) = 0.0;
            }
            
            if ( (m*l == 4) || (m*l == -4) ){
                HP.Element(2*i,2*j) += DE;
                HP.Element(2*i+1,2*j+1) += DE;
            }
        }
    }
    M.Product(HP,K);
    N.Product(HP,L);
    O.Product(M,M);
    P.Product(N,N);
    
    Gd_r.Product(K,O);
    Gd_a.Product(L,P);
    for(int i=0;i<Gd_r.Dim()*Gd_r.Dim();i++){
        Gd_r.Element(i) += K.Element(i);
        Gd_a.Element(i) += L.Element(i);
    }
#endif
}

double Hamiltonian::getNimp(){ return n_imp;}
complex<double> Hamiltonian::getGd(int m, int l, int s, int t, char control){
    if(control == 'r') return Gd_r.Element(2*m+s,2*l+t);
    else if (control == 'a') return Gd_a.Element(2*m+s,2*l+t);
    else return 0;
}


// Systems
MultiLayer_TBSQ::MultiLayer_TBSQ(){
    setLayerStructure();
    integral_type = 's';
}
void MultiLayer_TBSQ::setLayerStructure(){
    stringstream stm;
    stm << num_layer;
    
    //cout << "Ferromagnetic Metal / Insulator" << endl;
    structure.resize(LAYER_MAX);
    
    // FM/NI
    if ( H_TYPE == 00){
        for(int i=1;i<LAYER_MAX;i++)
            structure[i].soc_r=0.0;
        h_name = "FMNI";
    }
    else if ( H_TYPE == 04 ){
        for(int i=0;i<LAYER_MAX;i++){
            structure[i].t_ip = 0.0;
            if (i%2 == 1){
                    structure[i].soc_r = -1.0;
                    structure[i].ex_sd = -SUB_MAGNETIZATION;
                }
            }
        h_name = "2DDF";
    }
    h_name += stm.str();
    
}
void MultiLayer_TBSQ::getMatrix(Complex_SquareMatrix& A, double kx, double ky, double kz){
    A.SetMatrix(2*num_layer);
    
    for(int k=0;k<num_layer;k++){
        MonoLayer l = root * structure[k];
        SpinRotation R(ROTATION_AXIS,l.rad_p,l.rad_a);
        
        for(int i=0; i<2; i++){
            for(int j=0; j<2; j++){
                A.Element(2*k+i,2*k+j) = -2.0 * l.t_ip * (cos(kx)+cos(ky)) * Pauli[0][i][j]
                                        - l.ex_sd * R.getRot(i,j)
                                        + l.soc_r * (-sin(ky)*Pauli[1][i][j]+sin(kx)*Pauli[2][i][j]);
                if (k > 0){
                    A.Element(2*(k-1)+i,2*k+j) = -l.t_op*Pauli[0][i][j];
                    A.Element(2*k+i,2*(k-1)+j) = -l.t_op*Pauli[0][i][j];
                }
            }
        }
    }
}

void MultiLayer_TBSQ::getVelocity(Complex_SquareMatrix& A, double kx, double ky, double kz, char control){
    A.SetMatrix(2*num_layer);
    for(int k=0;k<num_layer;k++){
        MonoLayer l = root*structure[k];

        for(int i=0; i<2; i++){
            for(int j=0; j<2; j++){
                if(control == 'x')
                    A.Element(2*k+i,2*k+j) = 2 * l.t_ip * sin(kx) * Pauli[0][i][j] + l.soc_r * cos(kx)*Pauli[2][i][j];
                else
                    A.Element(2*k+i,2*k+j) = 2 * l.t_ip *sin(ky) * Pauli[0][i][j] + l.soc_r * (-cos(ky))*Pauli[1][i][j];
            }
        }
    }
}

void MultiLayer_TBSQ::getSpinVelocity(Complex_SquareMatrix& A, double kx, double ky, double kz, char control){
    A.SetMatrix(2*num_layer);
    for(int k=0;k<num_layer;k++){
        MonoLayer l = root*structure[k];

        for(int i=0; i<2; i++){
            for(int j=0; j<2; j++){
                if (control == 'x')
                    A.Element(2*k+i,2*k+j) = 2 * l.t_ip * sin(kx) * Pauli[3][i][j];
                else A.Element(2*k+i,2*k+j) = 2 * l.t_ip * sin(ky) * Pauli[3][i][j];
            }
        }
    }
}

void MultiLayer_TBSQ::getRashbaTorque(Complex_SquareMatrix& A, double kx, double ky, double kz, char control){
    A.SetMatrix(2*num_layer);
    for(int k=0;k<num_layer;k++){
        MonoLayer l = root*structure[k];

        for(int i=0; i<2; i++){
            for(int j=0; j<2; j++){
                if (control == 'x')
                    A.Element(2*k+i,2*k+j) = l.soc_r *  sin(kx)*Pauli[3][i][j];
                else if (control == 'y')
                    A.Element(2*k+i,2*k+j) = l.soc_r *  sin(ky)*Pauli[3][i][j];
                else 
                    A.Element(2*k+i,2*k+j) = -l.soc_r * ( sin(kx)*Pauli[1][i][j] + sin(ky)*Pauli[2][i][j] );
            }
        }
    }
}
 
void MultiLayer_TBSQ::getExchangeTorque(Complex_SquareMatrix& A, double kx, double ky, double kz, char control){
    A.SetMatrix(2*num_layer);
    for(int k=0;k<num_layer;k++){
        MonoLayer l = root*structure[k];

        for(int i=0; i<2; i++){
            for(int j=0; j<2; j++){
                if (control == 'x')
                    A.Element(2*k+i,2*k+j) = l.ex_sd*(-sin(l.rad_p)*sin(l.rad_a)*Pauli[3][i][j] + cos(l.rad_p)*Pauli[2][i][j]);
                else if (control == 'y')
                    A.Element(2*k+i,2*k+j) = l.ex_sd*(sin(l.rad_p)*cos(l.rad_a)*Pauli[3][i][j] - cos(l.rad_p)*Pauli[1][i][j]);
                else 
                    A.Element(2*k+i,2*k+j) = l.ex_sd*(sin(l.rad_p)*sin(l.rad_a)*Pauli[1][i][j] - sin(l.rad_p)*cos(l.rad_a)*Pauli[2][i][j]);
            }
        }
    }
}
void MultiLayer_TBSQ::getTransferTorque(Complex_SquareMatrix& A, double kx, double ky, double kz, char control){
    int spin;
    A.SetMatrix(2*num_layer);
    if (control == 'x') spin=1;
    else if (control == 'y') spin=2;
    else if (control == 'z') spin=3;
    else spin = 0;
    
    for(int k=0;k<num_layer;k++){
        MonoLayer l = root*structure[k];

        for(int i=0; i<2; i++){
            for(int j=0; j<2; j++){
                if (k != 0){
                    A.Element(2*(k-1)+i,2*k+j) = l.t_op/(2.0*I)*Pauli[spin][i][j];
                    A.Element(2*k+i,2*(k-1)+j) = -l.t_op/(2.0*I)*Pauli[spin][i][j];
                }
                if (k != num_layer-1){
                    A.Element(2*(k+1)+i,2*k+j) = l.t_op/(2.0*I)*Pauli[spin][i][j];
                    A.Element(2*k+i,2*(k+1)+j) = -l.t_op/(2.0*I)*Pauli[spin][i][j];
                }
            }
        }
    }
}

void MultiLayer_TBSQ::getSpinOperator(Complex_SquareMatrix& A, double kx, double ky, double kz, char control){
    A.SetMatrix(2*num_layer);
    for(int k=0;k<num_layer;k++){

        for(int i=0; i<2; i++){
            for(int j=0; j<2; j++){
                if (control == 'x')
                    A.Element(2*k+i,2*k+j) = Pauli[1][i][j];
                else if (control == 'y')
                    A.Element(2*k+i,2*k+j) = Pauli[2][i][j];
                else 
                    A.Element(2*k+i,2*k+j) = Pauli[3][i][j];
            }
        }
    }
}


void MultiLayer_TBSQ::getHybridization(vector<complex<double> > &A, double kx, double ky, double kz){
    if ( SD_SCATTERING != 0){
        cout << "TBSQ's Hds has not defined yet." << endl;
        exit(1);
    }
}

void MultiLayer_TBSQ::OutputStructure(ofstream& fout){
    fout << "# Tight-binding square lattice Hamiltonian (" << h_name << ": " << num_layer << " layer)"  << endl;
    fout << "# H[n]:, t_ip, EX_sd, RSOI, theta(dig), phi(dig), t_op" << endl;
    for(int i=0; i<num_layer; i++){
        int n=(num_layer-1)-i;
        MonoLayer l = root*structure[n];
        fout << "# H[" << n+1 << "], ";
        fout << l.t_ip << ", ";
        fout << l.ex_sd << ", ";
        fout << l.soc_r << ", ";
        fout << rad2dig(l.rad_p) << ", ";
        fout << rad2dig(l.rad_a);
        if (i != 0) fout << ", " << l.t_op;
        fout << endl;
    }
    fout << endl;
}


ElectronGas::ElectronGas(){
    //cout << "Electron gas" << endl;
    num_layer = 1;
    #if H_TYPE == 90
        integral_type = 'b';
        h_name = "3DEG";
        z_comp = 1.0;
    #else
        integral_type = 's';
        h_name = "2DEFG";
        z_comp = 0.0;
    #endif
}
void ElectronGas::getMatrix(Complex_SquareMatrix& A, double kx, double ky, double kz){
    A.SetMatrix(2);

    MonoLayer l = root;
    SpinRotation R(ROTATION_AXIS,l.rad_p,l.rad_a);
    
    for(int i=0; i<2; i++){
        for(int j=0; j<2; j++){
            A.Element(i,j) = l.t_ip * (kx*kx + ky*ky + z_comp*kz*kz) * Pauli[0][i][j]
                                    - l.ex_sd * R.getRot(i,j)
                                    + l.soc_r * (-ky*Pauli[1][i][j]+kx*Pauli[2][i][j]);
        }
    }
}

void ElectronGas::getVelocity(Complex_SquareMatrix& A, double kx, double ky, double kz, char control){
    A.SetMatrix(2);
    MonoLayer l = root;

    for(int i=0; i<2; i++){
        for(int j=0; j<2; j++){
            if(control == 'x')
                A.Element(i,j) = 2 * l.t_ip * kx * Pauli[0][i][j] + l.soc_r * Pauli[2][i][j];
            else if (control == 'y')
                A.Element(i,j) = 2 * l.t_ip * ky * Pauli[0][i][j] - l.soc_r * Pauli[1][i][j];
            else if (control == 'd')
                A.Element(i,j) = M_SQRT1_2*(2 * l.t_ip * (kx + ky) * Pauli[0][i][j] + l.soc_r * (Pauli[2][i][j]-Pauli[1][i][j]));
        }
    }
}

void ElectronGas::getSpinVelocity(Complex_SquareMatrix& A, double kx, double ky, double kz, char control){
    A.SetMatrix(2);
    MonoLayer l = root;

    for(int i=0; i<2; i++){
        for(int j=0; j<2; j++){
            if (control == 'x')
                A.Element(i,j) = 2 * l.t_ip * kx * Pauli[3][i][j];
            else
                A.Element(i,j) = 2 * l.t_ip * ky * Pauli[3][i][j];
        }
    }
}

void ElectronGas::getRashbaTorque(Complex_SquareMatrix& A, double kx, double ky, double kz, char control){
    A.SetMatrix(2);
    MonoLayer l = root;

    for(int i=0; i<2; i++){
        for(int j=0; j<2; j++){
            if (control == 'x')
                A.Element(i,j) = l.soc_r *  kx * Pauli[3][i][j];
            else if (control == 'y')
                A.Element(i,j) = l.soc_r *  ky * Pauli[3][i][j];
            else
                A.Element(i,j) = -l.soc_r * ( kx * Pauli[1][i][j] + ky * Pauli[2][i][j] );
        }
    }
}
 
void ElectronGas::getExchangeTorque(Complex_SquareMatrix& A, double kx, double ky, double kz, char control){
    A.SetMatrix(2);
    MonoLayer l = root;

    for(int i=0; i<2; i++){
        for(int j=0; j<2; j++){
            if (control == 'x')
                A.Element(i,j) = l.ex_sd*(-sin(l.rad_p)*sin(l.rad_a)*Pauli[3][i][j] + cos(l.rad_p)*Pauli[2][i][j]);
            else if (control == 'y')
                A.Element(i,j) = l.ex_sd*(sin(l.rad_p)*cos(l.rad_a)*Pauli[3][i][j] - cos(l.rad_p)*Pauli[1][i][j]);
            else
                A.Element(i,j) = l.ex_sd*(sin(l.rad_p)*sin(l.rad_a)*Pauli[1][i][j] - sin(l.rad_p)*cos(l.rad_a)*Pauli[2][i][j]);
        }
    }
}
void ElectronGas::getTransferTorque(Complex_SquareMatrix& A, double kx, double ky, double kz, char control){
    cout << "It hasn't defined yet." << endl;
    exit(1);
}

void ElectronGas::getSpinOperator(Complex_SquareMatrix& A, double kx, double ky, double kz, char control){
    A.SetMatrix(2);

    for(int i=0; i<2; i++){
        for(int j=0; j<2; j++){
            if (control == 'x')
                A.Element(i,j) = Pauli[1][i][j];
            else if (control == 'y')
                A.Element(i,j) = Pauli[2][i][j];
            else 
                A.Element(i,j) = Pauli[3][i][j];
        }
    }
}

void ElectronGas::getHybridization(vector<complex<double> > &A, double kx, double ky, double kz){
#if H_TYPE == 90
    for(int i=0; i<5; i++)
        A[i] = v_sd*SphericalHarmonicFunction(kx,ky,kz,i-2);
#elif H_TYPE == 91
    double phi_k = atan2(ky,kx);
    const double cY = sqrt(5.0/(6.0*M_PI));
    A[0] = cY *v_sd*exp(-2.0*I*phi_k);
    A[4] = cY *v_sd*exp(2.0*I*phi_k);
    A[1]=A[2]=A[3]=0.0;
#endif
}


void ElectronGas::OutputStructure(ofstream& fout){
    fout << "# Electron gas Hamiltonian (" << h_name <<  ")"  << endl;
    fout << "# H:, t, EX_sd, RSOI, theta(dig), phi(dig)" << endl;
    
    MonoLayer l = root;
    fout << "# H:, ";
    fout << l.t_ip << ", ";
    fout << l.ex_sd << ", ";
    fout << l.soc_r << ", ";
    fout << rad2dig(l.rad_p) << ", ";
    fout << rad2dig(l.rad_a);
    fout << endl;
    
    fout << endl;
}









































