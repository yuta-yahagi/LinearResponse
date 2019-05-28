// kubostreda.hpp

#ifndef KUBOSTREDA_HPP
#define KUBOSTREDA_HPP

#include "hamiltonian.hpp"


double surface_term(Complex_SquareMatrix& A, Complex_SquareMatrix& B, Complex_SquareMatrix& Gr, Complex_SquareMatrix& Ga);
double sea_term(Complex_SquareMatrix& A, Complex_SquareMatrix& B, Complex_SquareMatrix& Gr, Complex_SquareMatrix& Ga);
inline void GreenFunction(Hamiltonian& H, Complex_SquareMatrix& Hmat, double eps, double kx, double ky, double kz, char control);
inline void Born_approx(Hamiltonian& H, std::vector<std::complex<double> > &Vds, Complex_SquareMatrix& S, char control);

class Integral {
    public:
        double E_Integral(Hamiltonian& H, double eps);
        double K_Integral(Hamiltonian& H, double eps);
        double E_Integral(Hamiltonian& H, double eps, char control);
        double K_Integral(Hamiltonian& H, double eps, char control);
        virtual double Function(Hamiltonian& H, double eps, double kx, double ky, double kz, char control)=0;
};

class FermiLevel : public Integral{
    public:
        void calcFE(Hamiltonian& H);
        double getDeficiency();
        double getDOSatEF();
        double getEF();
        double Function(Hamiltonian& H, double eps, double kx, double ky, double kz, char control);

    private:
        double Energy_Iterator(Hamiltonian& H, double max, double min, double n_sum, int count);
        double EF;
        double DOSatEF;
        double deficiency;
        double n_max;
        double n_tail;
};

class DOS : public Integral {
    public:
        double Function(Hamiltonian& H, double eps, double kx, double ky, double kz, char control);
};

class DCC_Surface : public Integral {
    public:
        double Function(Hamiltonian& H, double eps, double kx, double ky, double kz, char control);
};

class AHC_Surface : public Integral {
    public:
        double Function(Hamiltonian& H, double eps, double kx, double ky, double kz, char control);
};

class AHC_Sea : public Integral {
    public:
        double Function(Hamiltonian& H, double eps, double kx, double ky, double kz, char control);
};

class SHC_Surface : public Integral {
    public:
        double Function(Hamiltonian& H, double eps, double kx, double ky, double kz, char control);
};

class SHC_Sea : public Integral {
    public:
        double Function(Hamiltonian& H, double eps, double kx, double ky, double kz, char control);
};

class RasT_Surface : public Integral {
    public:
        double Function(Hamiltonian& H, double eps, double kx, double ky, double kz, char control);
};
class RasT_Sea : public Integral {
    public:
        double Function(Hamiltonian& H, double eps, double kx, double ky, double kz, char control);
};

class ExcT_Surface : public Integral {
    public:
        double Function(Hamiltonian& H, double eps, double kx, double ky, double kz, char control);
};
class ExcT_Sea : public Integral {
    public:
        double Function(Hamiltonian& H, double eps, double kx, double ky, double kz, char control);
};
class STT_Surface : public Integral {
    public:
        double Function(Hamiltonian& H, double eps, double kx, double ky, double kz, char control);
};
class STT_Sea : public Integral {
    public:
        double Function(Hamiltonian& H, double eps, double kx, double ky, double kz, char control);
};
class Spin_Surface : public Integral {
    public:
        double Function(Hamiltonian& H, double eps, double kx, double ky, double kz, char control);
};
class Spin_Sea : public Integral {
    public:
        double Function(Hamiltonian& H, double eps, double kx, double ky, double kz, char control);
};

class Utot_Sea : public Integral {
    public:
        double Function(Hamiltonian& H, double eps, double kx, double ky, double kz, char control);
};

class Mz_Sea: public Integral {
    public:
        double Function(Hamiltonian& H, double eps, double kx, double ky, double kz, char control);
};

class Calculator {
    public:
        virtual void Run(Hamiltonian& H)=0;
        std::string Outputs();
        std::string Results();
        std::string Name();
        double Value();
        double Value(int i);
    protected:
        std::string name;
        std::vector<std::string> outputs;
        std::vector<double> results;
        double coefficient;
};

class Density_Of_States : public Calculator {
    public:
        Density_Of_States();
        void Run(Hamiltonian& H);
        double Imp_DOS(Hamiltonian &H, char control);
};

class DC_Conductivity : public Calculator {
    public:
        DC_Conductivity();
        void Run(Hamiltonian& H);
};

class Hall_Conductivity : public Calculator {
    public:
        Hall_Conductivity();
        void Run(Hamiltonian& H);
};

class SpinHall_Conductivity : public Calculator {
    public:
        SpinHall_Conductivity();
        void Run(Hamiltonian& H);
};

class SpinOrbitTorque : public Calculator {
    public:
        SpinOrbitTorque();
        void Run(Hamiltonian& H);
};

class SpinAccumulation : public Calculator {
    public:
        SpinAccumulation();
        void Run(Hamiltonian& H);
};

class Magnetization : public Calculator {
    public:
        Magnetization();
        void Run(Hamiltonian& H);
};

class TotalEnergy : public Calculator {
    public:
        TotalEnergy();
        void Run(Hamiltonian& H);
};

class MR_Amplitude : public Calculator {
    public:
        MR_Amplitude();
        void Run(Hamiltonian& H);
};

class Anisotropic_Energy : public Calculator {
    public:
        Anisotropic_Energy();
        void Run(Hamiltonian& H);
};

// Test tool
class BerryPhase : public Calculator {
    public:
        BerryPhase();
        void Run(Hamiltonian& H);
        double Energy(Hamiltonian& H, double kx, double ky, double kz, char sgn);
        double Fermi(double eps);
        double BzK(Hamiltonian& H, double kx, double ky, double kz,char sgn);
        void Error(double x);
};


#endif
