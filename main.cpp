#include<iostream>
#include<sstream>
#include<iomanip>
#include<cstdlib>
#include<time.h>
#include<omp.h>
#include "hamiltonian.hpp"
#include "config.hpp"
#include "kubostreda.hpp"


using namespace std;

string get_symbol(void);
string get_paramname(void);
string get_time(void);
string out_time(time_t t);

int main()
{
    clock_t start = clock();
    
    #if H_TYPE/10 < 5
        MultiLayer_TBSQ hamiltonian;
    #elif H_TYPE/10 < 9
        return -1;
    #elif H_TYPE/10 == 9
        ElectronGas hamiltonian;
    #else 
        return -1;
    #endif
    
    #if OUTPUT == 0
        DC_Conductivity calc;
    #elif OUTPUT == 1
        Hall_Conductivity calc;
    #elif OUTPUT == 2
        SpinHall_Conductivity calc;
    #elif OUTPUT == 3
        Density_Of_States calc;
    #elif OUTPUT == 4
        Magnetization calc;
    #elif OUTPUT == 5
        TotalEnergy calc;
    #elif OUTPUT == 6
        MR_Amplitude calc;
    #elif OUTPUT == 7
        SpinOrbitTorque calc;
    #elif OUTPUT == 8
        SpinAccumulation calc;
    #elif OUTPUT == 9
        Anisotropic_Energy calc;
    #else 
        return -1;
    #endif


    
    ofstream fout;
    string symb,iname;   
    symb = get_symbol();
    iname = get_paramname();

#if TEST_MODE == 0    
    stringstream fname,command,outfile,logfile;
    string home(getenv("HOME"));
    
    fname << PREFIX  << hamiltonian.H_Name() << "_" << calc.Name() << "vs" << symb  << SUFFIX << ".csv";
    command << "mkdir " << home << "/DATA/" << DIR_NAME;
    system(command.str().c_str());
    command.str("");
    command << "mkdir " << home << "/DATA/CALC_LOG/" << DIR_NAME;
    system(command.str().c_str());
    command.str("");
    
    logfile << home << "/DATA/CALC_LOG/" << DIR_NAME << "/" << get_time() << fname.str();
    outfile << home << "/DATA/" << DIR_NAME << "/" << fname.str();
    fout.open(logfile.str().c_str());
    if (!fout) {
        cout << "Cannot open \"" << logfile.str().c_str() << "\"." <<endl;
        return -1;
    }
    else cout << "Open file:\"" << logfile.str().c_str() << "\"." << endl;
    fout << "# Kubo-Streda formulation ver." << VERSION << " (Modified " << DATE << " Done " << get_time() << ")" << endl;
    fout << endl;
    hamiltonian.OutputHamiltonian(fout);
    fout << "# " << iname << calc.Outputs() << endl;
#elif TEST_MODE ==1
#warning "Test Mode"
    fout.open("result.csv");
    if (!fout) {
        cout << "Cannot open \"result.csv\"." <<endl;
        return -1;
    }
    cout << "# Kubo-Streda formulation (ver." << VERSION << ", modified " << DATE << ")" << endl;
#endif
    
    cout << "Calculation start." << endl;
    
    cout << "# " << iname << calc.Outputs() << endl;
    

    double min=INPUT_MIN, max=INPUT_MAX;
    double dx=(max-min)/(double)INPUT_DIV;
    for(int i=0; i<=INPUT_DIV; i++){
        double x = min + (double)i*dx;
        hamiltonian.setParam(x);
        //~ hamiltonian.printParam();
        
        calc.Run(hamiltonian);
        cout << x << ", ";
        cout << calc.Results();
        cout << endl;
        fout << x << ", ";
        fout << calc.Results();
        fout << endl;

    }
    
    fout.close();
#if TEST_MODE == 0
    command << "cp " << logfile.str() << " " << outfile.str();
    system(command.str().c_str());
#endif

	cout << endl << "Complete!" << endl;
    clock_t end = clock();
    cout << "Duration = " << out_time(end - start) << endl;
    
	return 0;
}
#if TEST_MODE == 2
int main(void){
    clock_t start = clock();
    
    cout << "test mode : " << TEST_MODE << endl;
    
    Hamiltonian h;
    Hall_Conductivity calc;
    BerryPhase prob;
    ofstream fout;
    fout.open("testlog.csv");
    if (!fout) {
        cout << "Cannot open file." <<endl;
        return -1;
    }
    fout << "# Accuracy test." << endl;
    double min=INPUT_MIN, max=INPUT_MAX;
    double dx=(max-min)/(double)INPUT_DIV;
    
    for(int i=0; i<=INPUT_DIV; i++){
        double x = min + (double)i*dx;
        h.setParam(x);
        calc.Run(h);
        prob.Run(h);
        prob.Error(calc.Value());
        cout << "AHC(" << x << "):" << calc.Value() << ", " << prob.Results() << endl;
        fout << x << ", "<< calc.Value() << ", " << prob.Results() << endl;
    }
	cout <<"Complete." << endl;
    fout.close();
    
    clock_t end = clock();
    cout << "Duration = " << out_time(end - start) << endl;
    
    return 0;
}

#elif TEST_MODE == 3
int main(void){
    clock_t start = clock();
    
    int sum=0;
    cout << "test mode : " << TEST_MODE << endl;
    
    #ifdef _OPENMP
        omp_set_num_threads(PARALLEL_NUM);
    
        ofstream fout;
        fout.open("omptest.csv");
        if (!fout) {
            cout << "Cannot open file." <<endl;
            return -1;
        }
    
        fout << "# PARALLEL_NUM = " << PARALLEL_NUM << endl;
        fout << "# process, thread_num" << endl;
        #pragma omp parallel for schedule(static) reduction (+:sum)
        for(int i=0; i<=200; i++){
            //~ cout << i << ", "<< omp_get_thread_num() << endl;
            #pragma omp critical
            fout << i << ", "<< omp_get_thread_num() << endl;
            sum += i;
        }
        cout << sum << endl;
        cout <<"Complete." << endl;
    #else
        # warning "Unlink OPEN_MP"
    #endif
    
    clock_t end = clock();
    cout << "Duration = " << out_time(end - start) << endl;
    
    return 0;
}
#elif TEST_MODE == 4
int main(void){
    clock_t start = clock();
    
    cout << "test mode : " << TEST_MODE << endl;
    
    Hamiltonian h;
    Complex_SquareMatrix Hmat;
    MR_Amplitude calc;
    double min=INPUT_MIN, max=INPUT_MAX;
    double dx=(max-min)/(double)INPUT_DIV;
    //~ for(int i=0; i<=INPUT_DIV; i++){
        //~ double x = min + (double)i*dx;
        //~ h.setParam(x);
        //~ calc.Run(h);
        h.OutputHamiltonian();
        h.getTransferTorque(Hmat,M_PI_2,0,'z');
        Hmat.PrintMatrix();
        //~ cout << x << ", " << calc.Value(5) <<  endl;
    //~ }
	cout <<"Complete." << endl;
    clock_t end = clock();
    cout << "Duration = " << out_time(end - start) << endl;
    
    return 0;
}
#endif

string get_symbol(void){
    switch (INPUT) {
        case 00:
            return "EXC";
        case 01:
            return "RSO";
        case 02:
            return "Tip";
        case 03:
            return "Top";
        case 04:
            if ( ROTATION_AXIS == 'z')
                return "PAz";
            else if ( ROTATION_AXIS == 'x')
                return "PAx";
            else return "PAy";
        case 05:
            if ( ROTATION_AXIS == 'z')
                return "AAz";
            else if ( ROTATION_AXIS == 'x')
                return "AAx";
            else return "AAy";
        case 06:
            return "sWID";
            
        case 10:
            return "HEX";
        case 11:
            return "ASO";
        case 12:
            return "Eimp";
        case 13:
            return "CFOh";
        case 14:
            return "dWID";
        case 15:
            return "CFTr";
            
        case 20:
            return "Nimp";
        case 21:
            return "Vsd";

        case 90:
            return "FEN";
        case 91:
            return "FIL";
        case 92:
            return "LYR";
    }
}
string get_paramname(void){
    switch (INPUT) {
        case 00:
            return "sd_Exchange";
        case 01:
            return "Rashba_SOI";
        case 02:
            return "Hopping_ip";
        case 03:
            return "Hopping_op";
        case 04:
            return "M_theta";
        case 05:
            return "M_phi";
        case 06:
            return "s_Lorentzian_width";
            
        case 10:
            return "Hund_coupling";
        case 11:
            return "Atomic_SOI";
        case 12:
            return "Impurity_level";
        case 13:
            return "Crystal_field(Oh)";
        case 14:
            return "d_Lorentzian_width";
        case 15:
            return "Tetragonal_crystal_field";
            
        case 20:
            return "Impurity_density";
        case 21:
            return "Hybridization";

        case 90:
            return "Fermi_level";
        case 91:
            return "Filling";
        case 92:
            return "Layer";
    }
}



string get_time(void)
{
    stringstream sstr;
    time_t     current;
    struct tm  *local;

    time(&current);                     // get current time
    local = localtime(&current);        // transform to localtime
    
    
    sstr << setw(2) << setfill('0') << local->tm_year-100;
    sstr << setw(2) << setfill('0') << local->tm_mon+1;
    sstr << setw(2) << setfill('0') << local->tm_mday << "_";
    sstr << setw(2) << setfill('0') << local->tm_hour;
    sstr << setw(2) << setfill('0') << local->tm_min;
    return sstr.str();
}

string out_time(time_t t){
    double sec = (double)t / CLOCKS_PER_SEC;
    double min = sec/60;
    double hour = min / 60;
    stringstream sstr;
    sstr << (int)hour << ":"  ;
    sstr << setw(2) << setfill('0') << (int)min % 60 << ":"  ;
    sstr << setw(2) << setfill('0') << (int)sec % 60;
    
    return sstr.str();
}












