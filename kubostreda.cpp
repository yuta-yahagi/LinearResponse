#include "kubostreda.hpp"
#include "config.hpp"

#include<iostream>
#include<sstream>
#include<iomanip>
#include<omp.h>
#include<float.h>
#include<stdlib.h>

using namespace std;



double surface_term(Complex_SquareMatrix& A, Complex_SquareMatrix& B, Complex_SquareMatrix& Gr, Complex_SquareMatrix& Ga){
    int dim = Gr.Dim();
    complex<double> sum(0.0,0.0);
    Complex_SquareMatrix K(dim),L(dim),M(dim),N(dim);
    
    // K = (Gr-Ga)
    K.Sub(Gr,Ga);
    // L = Gr*A
    L.Product(Gr,A);
    // M = B*Gr*A
    M.Product(B,L);
    // N = B*Gr*A*(Gr-Ga)
    N.Product(M,K);
    // sum = Tr[-B*Gr*A*(Gr-Ga)]
    sum -= N.Trace();
    
    // L = Ga*B
    L.Product(Ga,B);
    // M = A*Ga*B
    M.Product(A,L);
    // N = A*Ga*B*(Gr-Ga)
    N.Product(M,K);
    // sum = Tr[B*(Gr-Ga)*A*Ga - B*Gr*A*(Gr-Ga)]
    sum += N.Trace();
    
    return real(sum)/(4*M_PI);
}
double sea_term(Complex_SquareMatrix& A, Complex_SquareMatrix& B, Complex_SquareMatrix& Gr, Complex_SquareMatrix& Ga){
    int dim = Gr.Dim();
    complex<double> sum(0.0,0.0);
    Complex_SquareMatrix K(dim),L(dim),M(dim),N(dim);
    
    // K = A*Gr
    K.Product(A,Gr);
    // L = Gr*A*Gr
    L.Product(Gr,K);
    
    // M = B*Gr
    M.Product(B,Gr);
    // N = B*Gr*Gr*A*Gr
    N.Product(M,L);
    // sum += Tr[B*Gr*Gr*A*Gr]
    sum += N.Trace();
    
    // M = GrB
    M.Product(Gr,B);
    // N = Gr*B*Gr*A*Gr
    N.Product(M,L);
    // sum -= Tr[B*Gr*A*Gr*Gr]
    sum -= N.Trace();
    
    
    // K = A*Ga
    K.Product(A,Ga);
    // L = Ga*A*Ga
    L.Product(Ga,K);
    
    // M = B*Ga
    M.Product(B,Ga);
    // N = B*Ga*Ga*A*Ga
    N.Product(M,L);
    // sum -= Tr[B*Ga*Ga*A*Ga]
    sum -= N.Trace();

    // M = Ga*B
    M.Product(Ga,B);
    // N = Ga*B*Ga*A*Ga
    N.Product(M,L);
    // sum += Tr[B*Ga*A*Ga*Ga]
    sum += N.Trace();
    
    return real(sum)/(4*M_PI);
}

void GreenFunction(Hamiltonian& H, Complex_SquareMatrix& G, double eps, double kx, double ky, double kz, char control){
    H.getMatrix(G,kx,ky,kz);
    
    Complex_SquareMatrix SelfEnergy(G.Dim());
    std::complex<double> z;
    vector<complex<double> > Vds(5);
    
    H.getHybridization(Vds,kx,ky,kz);
    Born_approx(H,Vds,SelfEnergy,control);
    
    if (control == 'r') z = std::complex<double>(eps, H.getLifetime());
    else if (control == 'a') z = std::complex<double>(eps, -1*H.getLifetime());

    for (int i=0; i<G.Row(); i++)
        for (int j=0; j<G.Column(); j++)
            if (i==j) G.Element(i,i) = z-G.Element(i,i)-SelfEnergy.Element(i,j);
            else G.Element(i,j) = -G.Element(i,j)-SelfEnergy.Element(i,j);
    G.Inverse();
}

inline void Born_approx(Hamiltonian& H, vector<complex<double> > &Vds, Complex_SquareMatrix& S, char control){
    double n_imp = H.getNimp();
    for(int i=0;i<S.Row();i++)
        for(int j=0;j<S.Column();j++){
            S.Element(i,j)=0.0;
            for(int l=0;l<5;l++)
                for(int m=0;m<5;m++) 
                    S.Element(i,j) += n_imp*conj(Vds[m])*H.getGd(m,l,i,j,control)*Vds[l];
                    //~ if(i==j)S.Element(i,j) += conj(Vds[m])*Vds[l];
                    //~ if(i==j)S.Element(i,j) += imag(H.getGd(m,l,i,j,control));
        }
}


double Integral::E_Integral(Hamiltonian& H, double eps){ return E_Integral(H, eps, 0); }
double Integral::K_Integral(Hamiltonian& H, double eps){ return K_Integral(H, eps, 0); }

double Integral::E_Integral(Hamiltonian& H, double eps, char control){
    double sum=0.0;
    const double a=E_MIN, b=eps;
    const int div = E_DIV;
    const double de=(b-a)/(double)div;
    
    sum = (K_Integral(H,a) + K_Integral(H,b))*de*0.5;
    for(int i=1;i<div;i++){
        double x = a + i*de;

        sum += K_Integral(H,x,control)*de;
    }
    
    return sum;
}

double Integral::K_Integral(Hamiltonian& H, double eps, char control){
    double sum=0.0,coef;
    const double a=K_MIN, b=K_MAX;
    const int div=K_DIV;
    const double d=(b-a)/(double)div;
    char type;
    
    // Impurity
    H.setD_GreenFunc(eps);
    
    // Integral type
    type = H.check_Integral_type();
    if (type == 's')
        coef = 1.0/(4*M_PI*M_PI*(double)H.getLayer());
    else if (type == 'b')
        coef = 1.0/(8.0*M_PI*M_PI*M_PI);
    else {
        cout << "Integral type is undifined." << endl;
    }
        

    sum = 0.0;
#ifdef _OPENMP
  omp_set_num_threads(PARALLEL_NUM);
#pragma omp parallel for reduction(+:sum) schedule(static)
#endif
    for(int i=0; i<=div; i++){
        double ki = a + i*d;
        double dx = d;
        if (i == 0 || i == div) dx=0.5*d;
        
        for(int j=0; j<=div; j++){
            double kj = a + j*d;
            double dy = d;
            if (j == 0 || j == div) dy=0.5*d;
            
            
            if (type == 's')
                sum += Function(H, eps, ki, kj,0.0, control)*dx*dy;
            else if (type == 'b') {
                for (int l=0; l<=div; l++){
                    double kl = a + l*d;
                    double dz = d;
                    if (j == 0 || j == div) dz=0.5*d;
                    sum += Function(H, eps, ki, kj,kl, control)*dx*dy*dz;
                }
            }
            
            
            // ********* Test function *******
            //~ for(int l=0; l<=div; l++){
                //~ double kl = a + l*d;
                //~ double dz = d;
                //~ if (l == 0 || l == div) dz=0.5*d;
                //~ for(int m=0;m<5;m++)
                    //~ sum += dx*dy*dz;
            //~ }
        }
    }
    
    //~ sum = 0;
    //~ for(int m=0;m<5;m++)
        //~ if (control == 'u')sum += imag(H.getGd(m,m,0,0,'r'));
        //~ else sum += imag(H.getGd(m,m,1,1,'r'));
    //~ sum /= -M_PI*coef;
    
    return sum*coef;
}

void FermiLevel::calcFE(Hamiltonian& H){
    n_max = H.getFilling();
    if (INPUT == 90){
        EF = H.getFermiEnergy();
        DOSatEF = K_Integral(H, EF);
        deficiency = 0.0;
        n_tail = 0.0;
    }
    else if (EF_HALFWAY_MODE == 1){
        EF = INIT_FERMI_ENERGY;
        DOSatEF = K_Integral(H, EF);
        deficiency = 0.0;
        n_tail = 0.0;
    }
    else {
        n_tail = K_Integral(H,E_MIN)*(E_MAX-E_MIN)/(double)(E_DIV);
        EF = Energy_Iterator(H, E_MAX, E_MIN, 0.0, EF_ITERATION);
    }
}

double FermiLevel::Energy_Iterator(Hamiltonian& H, double max, double min, double n_sum, int count){
    double de = (max-min)/(double)E_DIV, eps,tmp;
    int i;
    
    count--;
    
    i=0;
    while (i <= E_DIV && n_sum < n_max){
        eps = min + i*de;
        tmp = K_Integral(H,eps);
        n_sum += tmp*de;
        i++;
    }
    if (i > E_DIV) {
        cout << "Setting error \"E_MAX\" : failed to explor Fermi level." << endl;
        exit(1);
    }
    
    n_sum -= tmp*de;
    if (count>0){
        return Energy_Iterator(H, eps+de, eps, n_sum, count);
    }
    else {
        DOSatEF = tmp;
        deficiency = n_max-n_sum;
        return eps;
    }
    
}

//~ void FermiLevel::calcFE(Hamiltonian& H){
    //~ double Nmax = H.getFilling(), Emax = E_MAX, Emin = E_MIN;
    //~ double dE = (Emax-Emin)/(double)E_DIV;
    //~ double Nsum = 0.0,eps,tmp,sgn;
    //~ 
    //~ if (INPUT == 90){
        //~ EF = H.getFermiEnergy();
        //~ DOSatEF = K_Integral(H, EF);
        //~ deficiency = 0.0;
    //~ }
    //~ else if (EF_HALFWAY_MODE == 1){
        //~ EF = INIT_FERMI_ENERGY;
        //~ DOSatEF = K_Integral(H, EF);
        //~ deficiency = 0.0;
    //~ }
    //~ else {
        //~ Nsum = K_Integral(H,Emin)*dE/2.0;
        //~ for (int i=1; i<=E_DIV; i++){
            //~ eps = Emin + dE*i;
        //~ 
            //~ tmp = K_Integral(H,eps);
            //~ if(Nsum+(tmp*dE) < Nmax)Nsum += tmp*dE;
            //~ else {
                //~ Nsum += (tmp*dE)/2.0;
                //~ break;
            //~ }
        //~ }
    //~ 
//~ // Re-iteration
        //~ if (Nsum > Nmax)sgn = -1.0;
        //~ else sgn = 1.0;
        //~ Emin = eps;
        //~ Emax = eps + sgn*dE;
        //~ dE = sgn*dE/(double)E_DIV;
    //~ 
        //~ Nsum += K_Integral(H,Emin)*dE/2.0;
        //~ for (int i=1; i<=E_DIV; i++){
            //~ eps = Emin + dE*i;
        //~ 
            //~ tmp = K_Integral(H,eps);
            //~ if(sgn*(Nsum+(tmp*dE)-Nmax) < 0)Nsum += tmp*dE;
            //~ else {
                //~ Nsum += (tmp*dE)/2.0;
                //~ break;
            //~ }
        //~ }
    //~ 
        //~ DOSatEF = tmp;
        //~ EF = eps;
        //~ deficiency = Nmax-Nsum;
    //~ }
//~ }


double FermiLevel::getEF(){ return EF; }
double FermiLevel::getDeficiency(){ return deficiency; }
double FermiLevel::getDOSatEF(){ return DOSatEF; }
double FermiLevel::Function(Hamiltonian& H, double eps, double kx, double ky, double kz, char control){
    Complex_SquareMatrix Gr;
    GreenFunction(H,Gr,eps,kx,ky,kz,'r');
    return -imag(Gr.Trace())/M_PI;
}

double DOS::Function(Hamiltonian& H, double eps, double kx, double ky, double kz, char control){
    double sum=0.0;
    int l=H.getLayer();
    Complex_SquareMatrix Gr;
    GreenFunction(H,Gr,eps,kx,ky,kz,'r');

    for(int i=0; i<l; i++)
        if (control == 'u')
            sum += imag(Gr.Element(2*i,2*i));
        else if (control == 'd')
            sum += imag(Gr.Element(2*i+1,2*i+1));
        else sum += (imag(Gr.Element(2*i,2*i))+imag(Gr.Element(2*i+1,2*i+1)));

    return -sum/M_PI;
}
double Utot_Sea::Function(Hamiltonian& H, double eps, double kx, double ky, double kz, char control){
    int l = H.getLayer();
    double sum=0.0;
    Complex_SquareMatrix Gr;
    H.getMatrix(Gr, kx, ky,kz);
    GreenFunction(H,Gr,eps,kx,ky,kz,'r');
    for (int i=0; i<l*2; i++)
        sum += -imag(Gr.Element(i,i))/M_PI;
    return sum*eps;
}

double Mz_Sea::Function(Hamiltonian& H, double eps, double kx, double ky, double kz, char control){
    int l = H.getLayer();
    double sum=0.0;
    Complex_SquareMatrix Gr;
    H.getMatrix(Gr, kx, ky, kz);
    GreenFunction(H,Gr,eps,kx,ky,kz,'r');
    for (int i=0; i<l; i++){
        sum += -imag(Gr.Element(2*i,2*i))/M_PI;
        sum -= -imag(Gr.Element(2*i+1,2*i+1))/M_PI;
    }
    return sum;
}

double DCC_Surface::Function(Hamiltonian& H, double eps, double kx, double ky, double kz, char control){
    Complex_SquareMatrix vx,Gr,Ga;
    if (NODAL_CONDUCTION == 1) H.getVelocity(vx,kx,ky,kz,'d');
    else H.getVelocity(vx,kx,ky,kz,'x');
    GreenFunction(H,Gr,eps,kx,ky,kz,'r');
    GreenFunction(H,Ga,eps,kx,ky,kz,'a');
    return surface_term(vx,vx,Gr,Ga);
}

double AHC_Surface::Function(Hamiltonian& H, double eps, double kx, double ky, double kz, char control){
    Complex_SquareMatrix vx,vy,Gr,Ga;
    H.getVelocity(vx,kx,ky,kz,'x');
    H.getVelocity(vy,kx,ky,kz,'y');
    GreenFunction(H,Gr,eps,kx,ky,kz,'r');
    GreenFunction(H,Ga,eps,kx,ky,kz,'a');
    return surface_term(vx,vy,Gr,Ga);
}
double AHC_Sea::Function(Hamiltonian& H, double eps, double kx, double ky, double kz, char control){
    Complex_SquareMatrix vx,vy,Gr,Ga;
    H.getVelocity(vx,kx,ky,kz,'x');
    H.getVelocity(vy,kx,ky,kz,'y');
    GreenFunction(H,Gr,eps,kx,ky,kz,'r');
    GreenFunction(H,Ga,eps,kx,ky,kz,'a');
    return sea_term(vx,vy,Gr,Ga);
}

double SHC_Surface::Function(Hamiltonian& H, double eps, double kx, double ky, double kz, char control){
    Complex_SquareMatrix A,B,Gr,Ga;
    H.getVelocity(A,kx,ky,kz,'x');
    H.getSpinVelocity(B,kx,ky,kz,'y');
    GreenFunction(H,Gr,eps,kx,ky,kz,'r');
    GreenFunction(H,Ga,eps,kx,ky,kz,'a');
    return surface_term(A,B,Gr,Ga);
}
double RasT_Surface::Function(Hamiltonian& H, double eps, double kx, double ky, double kz, char control){
    Complex_SquareMatrix A,B,Gr,Ga;
    H.getVelocity(A,kx,ky,kz,'x');
    H.getRashbaTorque(B,kx,ky,kz,control);
    GreenFunction(H,Gr,eps,kx,ky,kz,'r');
    GreenFunction(H,Ga,eps,kx,ky,kz,'a');
    return surface_term(A,B,Gr,Ga);
}
double RasT_Sea::Function(Hamiltonian& H, double eps, double kx, double ky, double kz, char control){
    Complex_SquareMatrix A,B,Gr,Ga;
    H.getVelocity(A,kx,ky,kz,'x');
    H.getRashbaTorque(B,kx,ky,kz,control);
    GreenFunction(H,Gr,eps,kx,ky,kz,'r');
    GreenFunction(H,Ga,eps,kx,ky,kz,'a');
    return sea_term(A,B,Gr,Ga);
}
double ExcT_Surface::Function(Hamiltonian& H, double eps, double kx, double ky, double kz, char control){
    Complex_SquareMatrix A,B,Gr,Ga;
    H.getVelocity(A,kx,ky,kz,'x');
    H.getExchangeTorque(B,kx,ky,kz,control);
    GreenFunction(H,Gr,eps,kx,ky,kz,'r');
    GreenFunction(H,Ga,eps,kx,ky,kz,'a');
    return surface_term(A,B,Gr,Ga);
}
double ExcT_Sea::Function(Hamiltonian& H, double eps, double kx, double ky, double kz, char control){
    Complex_SquareMatrix A,B,Gr,Ga;
    H.getVelocity(A,kx,ky,kz,'x');
    H.getExchangeTorque(B,kx,ky,kz,control);
    GreenFunction(H,Gr,eps,kx,ky,kz,'r');
    GreenFunction(H,Ga,eps,kx,ky,kz,'a');
    return sea_term(A,B,Gr,Ga);
}
double STT_Surface::Function(Hamiltonian& H, double eps, double kx, double ky, double kz, char control){
    Complex_SquareMatrix A,B,Gr,Ga;
    H.getVelocity(A,kx,ky,kz,'x');
    H.getTransferTorque(B,kx,ky,kz,control);
    GreenFunction(H,Gr,eps,kx,ky,kz,'r');
    GreenFunction(H,Ga,eps,kx,ky,kz,'a');
    return surface_term(A,B,Gr,Ga);
}
double STT_Sea::Function(Hamiltonian& H, double eps, double kx, double ky, double kz, char control){
    Complex_SquareMatrix A,B,Gr,Ga;
    H.getVelocity(A,kx,ky,kz,'x');
    H.getTransferTorque(B,kx,ky,kz,control);
    GreenFunction(H,Gr,eps,kx,ky,kz,'r');
    GreenFunction(H,Ga,eps,kx,ky,kz,'a');
    return sea_term(A,B,Gr,Ga);
}
double Spin_Surface::Function(Hamiltonian& H, double eps, double kx, double ky, double kz, char control){
    Complex_SquareMatrix A,B,Gr,Ga;
    H.getVelocity(A,kx,ky,kz,'x');
    H.getSpinOperator(B,kx,ky,kz,control);
    GreenFunction(H,Gr,eps,kx,ky,kz,'r');
    GreenFunction(H,Ga,eps,kx,ky,kz,'a');
    return surface_term(A,B,Gr,Ga);
}
double Spin_Sea::Function(Hamiltonian& H, double eps, double kx, double ky, double kz, char control){
    Complex_SquareMatrix A,B,Gr,Ga;
    H.getVelocity(A,kx,ky,kz,'x');
    H.getSpinOperator(B,kx,ky,kz,control);
    GreenFunction(H,Gr,eps,kx,ky,kz,'r');
    GreenFunction(H,Ga,eps,kx,ky,kz,'a');
    return sea_term(A,B,Gr,Ga);
}

// Deffinition of Calculator
string Calculator::Name(){return name; }
string Calculator::Outputs(){
    stringstream sstr;
    for(unsigned int i=0; i<outputs.size(); i++) sstr << ", " << outputs[i];
    return sstr.str();
}
string Calculator::Results(){
    stringstream sstr;
    sstr << scientific << setprecision(PRECISION) << uppercase;
    for(unsigned int i=0; i<results.size(); i++) sstr << results[i] << ", ";
    return sstr.str();
}
double Calculator::Value(){ return results[0]; }
double Calculator::Value(int i){ return results[i]; }


// DOS
Density_Of_States::Density_Of_States(){
    coefficient = 1.0;
    outputs.push_back("DOS [/t]");
    outputs.push_back("UP_DOS");
    outputs.push_back("DN_DOS");
    outputs.push_back("Imp_DOS [/t]");
    outputs.push_back("Imp_UP_DOS");
    outputs.push_back("Imp_DN_DOS");
    if (INPUT == 90){
        outputs.push_back("n_tot");
        outputs.push_back("n_up");
        outputs.push_back("n_dn");
    }
    else outputs.push_back("EF");
    name = "DOS";
}
void Density_Of_States::Run(Hamiltonian& H){
    FermiLevel EF_obj;
    DOS dos_obj;
    double fl,upd,dnd,impup,impdn;
    std::vector<double> temp;
    static double n_up,n_dn;
    const double de = (INPUT_MAX-INPUT_MIN)/(double)INPUT_DIV;
    
    EF_obj.calcFE(H);
    fl = EF_obj.getEF();
    H.setFermiEnergy(fl);
    results.push_back(EF_obj.getDOSatEF());
    
    upd = dos_obj.K_Integral(H, fl, 'u');
    dnd = dos_obj.K_Integral(H, fl, 'd');
    temp.push_back(upd+dnd);
    temp.push_back(upd);
    temp.push_back(dnd);
    
    impup = Imp_DOS(H,'u');
    impdn = Imp_DOS(H,'d');
    temp.push_back(impup+impdn);
    temp.push_back(impup);
    temp.push_back(impdn);
    if (INPUT == 90) {
        n_up += upd*de;
        n_dn += dnd*de;
        temp.push_back(n_up+n_dn);
        temp.push_back(n_up);
        temp.push_back(n_dn);
    }
    else temp.push_back(fl);
    results = temp;
}
double Density_Of_States::Imp_DOS(Hamiltonian &H, char control){
    double sum = 0;
    for(int m=0;m<5;m++)
        if (control == 'u')sum += imag(H.getGd(m,m,0,0,'r'));
        else sum += imag(H.getGd(m,m,1,1,'r'));
    sum /= -M_PI;
    return sum;
}

// DC conductivity
DC_Conductivity::DC_Conductivity(){
    coefficient = 2.0*M_PI;
    if ( NODAL_CONDUCTION ==1 ) outputs.push_back("sigma_dd [e^2/h]");
    else outputs.push_back("sigma_xx [e^2/h]");
    outputs.push_back("Resistivity [h/e^2]");
    outputs.push_back("EF");
    outputs.push_back("DOS(EF)");
    outputs.push_back("Deficiency of N");
    if ( NODAL_CONDUCTION ==1 ) name = "NODAL_DCC";
    else name = "DCC";
}
void DC_Conductivity::Run(Hamiltonian& H){
    FermiLevel EF_obj;
    DCC_Surface Sgm_obj;
    double fl,dos,def,s;
    std::vector<double> temp;
    
    EF_obj.calcFE(H);
    fl = EF_obj.getEF();
    H.setFermiEnergy(fl);
    dos=EF_obj.getDOSatEF();
    def = EF_obj.getDeficiency();
    
    s = coefficient*Sgm_obj.K_Integral(H,fl);
    
    temp.push_back(s);
    temp.push_back(1.0/s);
    temp.push_back(fl);
    temp.push_back(dos);
    temp.push_back(def);
    
    results = temp;
}

// Anomalous Hall Effect
Hall_Conductivity::Hall_Conductivity(){
    coefficient = 2.0*M_PI;
    outputs.push_back("Hall_conductivity [e^2/h]");
    outputs.push_back("surface_term");
    outputs.push_back("sea_term");
    outputs.push_back("EF");
    outputs.push_back("DOS(EF)");
    outputs.push_back("Deficiency of N");
    name = "AHC";
}
void Hall_Conductivity::Run(Hamiltonian& H){
    FermiLevel EF_obj;
    AHC_Surface Sgm1_obj;
    AHC_Sea Sgm2_obj;
    double fl,dos,def,s1=0.0,s2=0.0;
    std::vector<double> temp;
    
    EF_obj.calcFE(H);
    fl = EF_obj.getEF();
    H.setFermiEnergy(fl);
    dos=EF_obj.getDOSatEF();
    def = EF_obj.getDeficiency();
    
    s1 = Sgm1_obj.K_Integral(H,fl);
    if(SF_HALFWAY_MODE==1)
        s2=0;
    else
        s2 = Sgm2_obj.E_Integral(H,fl);
    s1*=coefficient;
    s2*=coefficient;
    
    temp.push_back(s1+s2);   
    temp.push_back(s1);
    temp.push_back(s2);
    temp.push_back(fl);
    temp.push_back(dos);
    temp.push_back(def);
    
    results=temp;
}

// Spin-Hall
SpinHall_Conductivity::SpinHall_Conductivity(){
    coefficient = -1.0/2.0;
    outputs.push_back("SpinHall_conductivity [e]");
    outputs.push_back("EF");
    outputs.push_back("DOS(EF)");
    outputs.push_back("Deficiency of N");
    name = "SHC";
}
void SpinHall_Conductivity::Run(Hamiltonian& H){
    FermiLevel EF_obj;
    SHC_Surface Sgm1_obj;
    double fl,dos,def,s1=0.0;
    std::vector<double> temp;
    
    EF_obj.calcFE(H);
    fl = EF_obj.getEF();
    H.setFermiEnergy(fl);
    dos=EF_obj.getDOSatEF();
    def = EF_obj.getDeficiency();
    
    s1 = coefficient*Sgm1_obj.K_Integral(H,fl);
    
    temp.push_back(s1);
    temp.push_back(fl);
    temp.push_back(dos);
    temp.push_back(def);
    
    results=temp;
}

// Spin-Orbit Torque
SpinOrbitTorque::SpinOrbitTorque(){
    coefficient = -1.0/2.0;
    if ( (MU_DIR==0) || (MU_DIR==1)){
        outputs.push_back("Torque_x [eE]");
        outputs.push_back("Rashba_Tx");
        outputs.push_back("Excange_Tx");
        outputs.push_back("STT_x");
    }
    if ((MU_DIR==0) || (MU_DIR==2)){
        outputs.push_back("Torque_y [eE]");
        outputs.push_back("Rashba_Ty");
        outputs.push_back("Excange_Ty");
        outputs.push_back("STT_y");
    }
    if ((MU_DIR==0) || (MU_DIR==3)){
        outputs.push_back("Torque_z [eE]");
        outputs.push_back("Rashba_Tz");
        outputs.push_back("Excange_Tz");
        outputs.push_back("STT_z");
    }
    outputs.push_back("EF");
    outputs.push_back("DOS(EF)");
    outputs.push_back("Deficiency of N");
    name = "SOT";
}
void SpinOrbitTorque::Run(Hamiltonian& H){
    FermiLevel EF_obj;
    RasT_Surface RST1_obj;
    RasT_Sea RST2_obj;
    ExcT_Surface EXT1_obj;
    ExcT_Sea EXT2_obj;
    STT_Surface STT1_obj;
    STT_Sea STT2_obj;
    
    int layer = H.getLayer();
    double fl,dos,def;
    std::vector<double> temp;
    
    EF_obj.calcFE(H);
    fl = EF_obj.getEF();
    H.setFermiEnergy(fl);
    dos=EF_obj.getDOSatEF();
    def = EF_obj.getDeficiency();
    
    if ( (MU_DIR==0) || (MU_DIR==1)){
        double s1=0.0,s2=0.0,s3=0.0;
        s1 = RST1_obj.K_Integral(H,fl,'x');
        s2 = EXT1_obj.K_Integral(H,fl,'x');
        if (layer > 1) s3 = STT1_obj.K_Integral(H,fl,'x');
        if(SF_HALFWAY_MODE==0){
            s1 += RST2_obj.E_Integral(H,fl,'x');
            s2 += EXT2_obj.E_Integral(H,fl,'x');
            if (layer > 1) s3 += STT2_obj.E_Integral(H,fl,'x');
        }
        s1*=coefficient;
        s2*=coefficient;
        s3*=coefficient;
        temp.push_back(s1+s2+s3);
        temp.push_back(s1);
        temp.push_back(s2);
        temp.push_back(s3);
    }
    if ((MU_DIR==0) || (MU_DIR==2)){
        double s1=0.0,s2=0.0,s3=0.0;
        s1 = RST1_obj.K_Integral(H,fl,'y');
        s2 = EXT1_obj.K_Integral(H,fl,'y');
        if (layer > 1) s3 = STT1_obj.K_Integral(H,fl,'y');
        if(SF_HALFWAY_MODE==0){
            s1 += RST2_obj.E_Integral(H,fl,'y');
            s2 += EXT2_obj.E_Integral(H,fl,'y');
            if (layer > 1) s3 += STT2_obj.E_Integral(H,fl,'y');
        }
        s1*=coefficient;
        s2*=coefficient;
        s3*=coefficient;
        temp.push_back(s1+s2+s3);
        temp.push_back(s1);
        temp.push_back(s2);
        temp.push_back(s3);
    }
    if ((MU_DIR==0) || (MU_DIR==3)){
        double s1=0.0,s2=0.0,s3=0.0;
        s1 = RST1_obj.K_Integral(H,fl,'z');
        s2 = EXT1_obj.K_Integral(H,fl,'z');
        if (layer > 1) s3 = STT1_obj.K_Integral(H,fl,'z');
        if(SF_HALFWAY_MODE==0){
            s1 += RST2_obj.E_Integral(H,fl,'z');
            s2 += EXT2_obj.E_Integral(H,fl,'z');
            if (layer > 1) s3 += STT2_obj.E_Integral(H,fl,'z');
        }
        s1*=coefficient;
        s2*=coefficient;
        s3*=coefficient;
        temp.push_back(s1+s2+s3);
        temp.push_back(s1);
        temp.push_back(s2);
        temp.push_back(s3);
    }
    temp.push_back(fl);
    temp.push_back(dos);
    temp.push_back(def);
    
    results = temp;
}

// Spin Accumulation
SpinAccumulation::SpinAccumulation(){
    coefficient = -1.0/(4.0*M_PI);
    if ((MU_DIR==0) || (MU_DIR==1)) outputs.push_back("S_x [heE]");
    if ((MU_DIR==0) || (MU_DIR==2)) outputs.push_back("S_y [heE]");
    if ((MU_DIR==0) || (MU_DIR==3)) outputs.push_back("S_z [heE]");
    outputs.push_back("EF");
    outputs.push_back("DOS(EF)");
    outputs.push_back("Deficiency of N");
    name = "REE";
}
void SpinAccumulation::Run(Hamiltonian& H){
    FermiLevel EF_obj;
    Spin_Surface S1_obj;
    Spin_Sea S2_obj;
    double fl,dos,def;
    std::vector<double> temp;
    
    EF_obj.calcFE(H);
    fl = EF_obj.getEF();
    H.setFermiEnergy(fl);
    dos=EF_obj.getDOSatEF();
    def = EF_obj.getDeficiency();
    
    if ( (MU_DIR==0) || (MU_DIR==1)){
        double s=0.0;
        s = coefficient*S1_obj.K_Integral(H,fl,'x');
        temp.push_back(s);
    }
    if ((MU_DIR==0) || (MU_DIR==2)){
        double s=0.0;
        s = coefficient*S1_obj.K_Integral(H,fl,'y');
        if(SF_HALFWAY_MODE==0) s += S2_obj.E_Integral(H,fl,'y');
        temp.push_back(s);
    }
    if ((MU_DIR==0) || (MU_DIR==3)){
        double s=0.0;
        s = coefficient*S1_obj.K_Integral(H,fl,'z');
        temp.push_back(s);
    }
    temp.push_back(fl);
    temp.push_back(dos);
    temp.push_back(def);
    
    results = temp;
}

// Magnetization
Magnetization::Magnetization(){
    coefficient = 1.0;
    outputs.push_back("Magnetization Mz [mu_B]");
    outputs.push_back("EF");
    outputs.push_back("DOS(EF)");
    outputs.push_back("Deficiency of N");
    name = "MAG";
}
void Magnetization::Run(Hamiltonian& H){
    FermiLevel EF_obj;
    Mz_Sea Mz_obj;
    double fl,dos,def,Mz=0.0;
    std::vector<double> temp;
    
    EF_obj.calcFE(H);
    fl = EF_obj.getEF();
    H.setFermiEnergy(fl);
    dos=EF_obj.getDOSatEF();
    def = EF_obj.getDeficiency();

    Mz = Mz_obj.E_Integral(H,fl);
    temp.push_back(Mz);
    temp.push_back(fl);
    temp.push_back(dos);
    temp.push_back(def);
    
    results = temp;
}


// Total Energy
TotalEnergy::TotalEnergy(){
    coefficient = 1.0;
    outputs.push_back("Total Energy Utot/t [.]");
    outputs.push_back("EF");
    outputs.push_back("DOS(EF)");
    outputs.push_back("Deficiency of N");
    name = "ENG";
}
void TotalEnergy::Run(Hamiltonian& H){
    FermiLevel EF_obj;
    Utot_Sea U_obj;
    double fl,dos,def,U=0.0;
    std::vector<double> temp;

    EF_obj.calcFE(H);
    fl = EF_obj.getEF();
    H.setFermiEnergy(fl);
    dos=EF_obj.getDOSatEF();
    def = EF_obj.getDeficiency();
    
    U = U_obj.E_Integral(H,fl);
    temp.push_back(U);
    temp.push_back(fl);
    temp.push_back(dos);
    temp.push_back(def);
    
    results = temp;
}

MR_Amplitude::MR_Amplitude(){
    coefficient = 2.0*M_PI;
    
    outputs.push_back("rho_x[R_K]");
    outputs.push_back("rho_y");
    outputs.push_back("rho_z");
    if(FOURFOLD==1){
        outputs.push_back("rho_xy(fourfold)");
        outputs.push_back("rho_yz");
        outputs.push_back("rho_zx");
        
        outputs.push_back("r0_xy");
        outputs.push_back("r0_yz");
        outputs.push_back("r0_zx");
        
        outputs.push_back("r2_xy");
        outputs.push_back("r2_yz");
        outputs.push_back("r2_zx");
        
        outputs.push_back("r4_xy");
        outputs.push_back("r4_yz");
        outputs.push_back("r4_zx");
    }
    
    else {
        outputs.push_back("rho_0");
        outputs.push_back("MRxy [%]");
        outputs.push_back("MRzy");
        outputs.push_back("MRzx");
    }
    name = "AMR";
}
void MR_Amplitude::Run(Hamiltonian& H){
    FermiLevel EF_obj;
    DCC_Surface Sgm_obj;
    double fl,s,Rx,Ry,Rz,R0;
    std::vector<double> temp;
    
    H.setDirection('x');
    EF_obj.calcFE(H);
    fl = EF_obj.getEF();
    H.setFermiEnergy(fl);
    s = coefficient*Sgm_obj.K_Integral(H,fl);
    Rx = 1.0/s;
    
    H.setDirection('y');
    EF_obj.calcFE(H);
    fl = EF_obj.getEF();
    H.setFermiEnergy(fl);
    s = coefficient*Sgm_obj.K_Integral(H,fl);
    Ry = 1.0/s;

    H.setDirection('z');
    EF_obj.calcFE(H);
    fl = EF_obj.getEF();
    H.setFermiEnergy(fl);
    s = coefficient*Sgm_obj.K_Integral(H,fl);
    Rz = 1.0/s;
    
    if(FOURFOLD==1){
        double Rxy,Ryz,Rzx;
        double rxy[3],ryz[3],rzx[3];
        
        H.setDirection(90,45);
        EF_obj.calcFE(H);
        fl = EF_obj.getEF();
        H.setFermiEnergy(fl);
        s = coefficient*Sgm_obj.K_Integral(H,fl);
        Rxy = 1.0/s;
        
        H.setDirection(45,90);
        EF_obj.calcFE(H);
        fl = EF_obj.getEF();
        H.setFermiEnergy(fl);
        s = coefficient*Sgm_obj.K_Integral(H,fl);
        Ryz = 1.0/s;
        
        H.setDirection(45,0);
        EF_obj.calcFE(H);
        fl = EF_obj.getEF();
        H.setFermiEnergy(fl);
        s = coefficient*Sgm_obj.K_Integral(H,fl);
        Rzx = 1.0/s;
        
        // Symmetrically decompose
        rxy[0] = 0.25*(Rx+Ry+2*Rxy);
        rxy[1] = 0.5*(Rx-Ry);
        rxy[2] = 0.25*(Rx+Ry-2*Rxy);
        
        ryz[0] = 0.25*(Ry+Rz+2*Ryz);
        ryz[1] = 0.5*(Ry-Rz);
        ryz[2] = 0.25*(Ry+Rz-2*Ryz);
        
        rzx[0] = 0.25*(Rz+Rx+2*Rzx);
        rzx[1] = 0.5*(Rz-Rx);
        rzx[2] = 0.25*(Rz+Rx-2*Rzx);
        
        temp.push_back(Rx);
        temp.push_back(Ry);
        temp.push_back(Rz);
        temp.push_back(Rxy);
        temp.push_back(Ryz);
        temp.push_back(Rzx);
        
        temp.push_back(rxy[0]);
        temp.push_back(ryz[0]);
        temp.push_back(rzx[0]);

        temp.push_back(rxy[1]);
        temp.push_back(ryz[1]);
        temp.push_back(rzx[1]);
        
        temp.push_back(rxy[2]);
        temp.push_back(ryz[2]);
        temp.push_back(rzx[2]);
    }
    else {
        R0=(Rx+Ry+Rz)/3.0;
        
        temp.push_back(Rx);
        temp.push_back(Ry);
        temp.push_back(Rz);
        temp.push_back(R0);
        temp.push_back((Rx-Ry)/R0*100);
        temp.push_back((Rz-Ry)/R0*100);
        temp.push_back((Rz-Rx)/R0*100);
    }
    results = temp;
}

Anisotropic_Energy::Anisotropic_Energy(){
    
    coefficient = 1.0;
    outputs.push_back("U0 [t]");
    outputs.push_back("Ux");
    outputs.push_back("Uy");
    outputs.push_back("Uz");
    outputs.push_back("Ux-Uy");
    outputs.push_back("Uz-Uy");
    outputs.push_back("Uz-Ux");
    name = "AEN";
}
void Anisotropic_Energy::Run(Hamiltonian& H){
    FermiLevel EF_obj;
    Utot_Sea U_obj;
    double fl,U,Ux,Uy,Uz;
    std::vector<double> temp;

    H.setDirection('x');
    EF_obj.calcFE(H);
    fl = EF_obj.getEF();
    H.setFermiEnergy(fl);
    Ux = U_obj.E_Integral(H,fl);
    
    H.setDirection('y');
    EF_obj.calcFE(H);
    fl = EF_obj.getEF();
    H.setFermiEnergy(fl);
    Ux = U_obj.E_Integral(H,fl);

    H.setDirection('z');
    EF_obj.calcFE(H);
    fl = EF_obj.getEF();
    H.setFermiEnergy(fl);
    Ux = U_obj.E_Integral(H,fl);
    
    U=(Ux+Uy+Uz)/3.0;
    
    temp.push_back(U);
    temp.push_back(Ux);
    temp.push_back(Uy);
    temp.push_back(Uz);
    temp.push_back(Ux-Uy);
    temp.push_back(Uz-Uy);
    temp.push_back(Uz-Ux);
    
    results = temp;
}



