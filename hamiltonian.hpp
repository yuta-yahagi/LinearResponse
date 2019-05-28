// hamiltonian.hpp
#ifndef HAMILTONIAN_HPP
#define HAMILTONIAN_HPP

#include"matrix.hpp"
#include<string>
#include<fstream>

const std::complex<double> I(0.0,1.0);
const std::complex<double> Pauli[4][2][2]={
  {{1.0+I*0.0,0.0+I*0.0},
   {0.0+I*0.0,1.0+I*0.0}},

  {{0.0+I*0.0,1.0+I*0.0},
   {1.0+I*0.0,0.0+I*0.0}},

  {{0.0+I*0.0,0.0-I*1.0},
   {0.0+I*1.0,0.0+I*0.0}},

  {{1.0+I*0.0,0.0+I*0.0},
   {0.0+I*0.0,-1.0+I*0.0}}
};
double dig2rad(double dig);
double rad2dig(double dig);
std::complex<double> SphericalHarmonicFunction(double x,double y,double z,int m);


struct MonoLayer {
    public:
        MonoLayer();
        double t_ip;
        double ex_sd;
        double soc_r;
        double rad_p;
        double rad_a;
        double t_op;
        const MonoLayer operator *(const MonoLayer& rhs) const;
        MonoLayer& operator*=(const MonoLayer& rhs);
};
struct Impurity {
    public:
        Impurity();
        double d_orbital_level;
        double d_orbital_width;
        double d_orbital_exchange;
        double d_orbital_soi;
        double d_orbital_cf;
        double d_orbital_tetra;
        double d_orbital_screened_charge;
};

class SpinRotation {
    public:
        SpinRotation();
        SpinRotation(double theta, double phi);
        SpinRotation(char axis, double theta, double phi);
        
        std::complex<double> getRot(int i, int j);
        
    private:
        std::complex<double> rotmat[2][2];
};

class Hamiltonian {
    public:
        Hamiltonian(); // FM/NI, FM/NM, etc...
        void printParam();
        void printParam(char option);
        char setParam(double a);
        void setFermiEnergy(double x);
        void setDirection(char dir);
        void setDirection(double theta, double phi);
        int getLayer();
        MonoLayer getRootParam();
        double getLifetime();
        double getFermiEnergy();
        double getFilling();
        char check_Integral_type(); // 3D: b, 2D: s, QW: w
        
        virtual void getMatrix(Complex_SquareMatrix& A, double kx, double ky, double kz)=0;
        virtual void getVelocity(Complex_SquareMatrix& A, double kx, double ky, double kz, char control)=0;
        virtual void getSpinVelocity(Complex_SquareMatrix& A, double kx, double ky, double kz,char control)=0;
        virtual void getRashbaTorque(Complex_SquareMatrix& A, double kx, double ky, double kz, char control)=0;
        virtual void getExchangeTorque(Complex_SquareMatrix& A, double kx, double ky, double kz, char control)=0;
        virtual void getTransferTorque(Complex_SquareMatrix& A, double kx, double ky, double kz, char control)=0;
        virtual void getSpinOperator(Complex_SquareMatrix& A, double kx, double ky, double kz, char control)=0;
        virtual void getHybridization(std::vector<std::complex<double> > &A,double kx, double ky, double kz)=0;
        
        // Impurity
        void setSD_Impurity();
        void setD_GreenFunc(double eps);
        double getNimp();
        std::complex<double> getGd(int m, int l, int s, int t, char control);
        
        std::string H_Name();
        //void OutputHamiltonian();
        void OutputHamiltonian(std::ofstream& fout);
        
    protected:
        void Initialize();
        int num_layer;
        std::string h_name;
        MonoLayer root;
        std::vector<MonoLayer> structure;
        double gamma;
        double epsf;
        double filling;
        double n_imp;
        char integral_type;
        
        double v_sd;
        
        // Impurity
        Impurity imp;
        Complex_SquareMatrix Hd;
        Complex_SquareMatrix Gd_r;
        Complex_SquareMatrix Gd_a;

        virtual void OutputStructure(std::ofstream& fout) = 0;
};

class MultiLayer_TBSQ : public Hamiltonian {
    public:
        MultiLayer_TBSQ();
        void getMatrix(Complex_SquareMatrix& A, double kx, double ky, double kz);
        void getVelocity(Complex_SquareMatrix& A, double kx, double ky, double kz, char control);
        void getSpinVelocity(Complex_SquareMatrix& A, double kx, double ky, double kz,char control);
        void getRashbaTorque(Complex_SquareMatrix& A, double kx, double ky, double kz, char control);
        void getExchangeTorque(Complex_SquareMatrix& A, double kx, double ky, double kz, char control);
        void getTransferTorque(Complex_SquareMatrix& A, double kx, double ky, double kz, char control);
        void getSpinOperator(Complex_SquareMatrix& A, double kx, double ky, double kz, char control);
        
        void getHybridization(std::vector<std::complex<double> > &A, double kx, double ky, double kz);
    
    private:
        void setLayerStructure();
        void OutputStructure(std::ofstream& fout);
};

class ElectronGas : public Hamiltonian {
    public:
        ElectronGas();
        void getMatrix(Complex_SquareMatrix& A, double kx, double ky, double kz);
        void getVelocity(Complex_SquareMatrix& A, double kx, double ky, double kz, char control);
        void getSpinVelocity(Complex_SquareMatrix& A, double kx, double ky, double kz,char control);
        void getRashbaTorque(Complex_SquareMatrix& A, double kx, double ky, double kz, char control);
        void getExchangeTorque(Complex_SquareMatrix& A, double kx, double ky, double kz, char control);
        void getTransferTorque(Complex_SquareMatrix& A, double kx, double ky, double kz, char control);
        void getSpinOperator(Complex_SquareMatrix& A, double kx, double ky, double kz, char control);
        
        void getHybridization(std::vector<std::complex<double> > &A, double kx, double ky, double kz);
    
    private:
        double z_comp;
        void OutputStructure(std::ofstream& fout);
};

#endif
