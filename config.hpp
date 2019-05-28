// config.hpp
#ifndef CONFIG_HPP
#define CONFIG_HPP

// Last modified date
#define DATE "2018/12/19"
#define VERSION "2.8"

// Name of output directory
#define DIR_NAME "DEFAULT"

// Prefix & suffix name of output file
#define PREFIX ""
#define SUFFIX "_Tetra"

// Number of thread for parallel
#define PARALLEL_NUM 4

// Control main func
// 0: Off
// 1: Console & Local Directory Output Mode
// 2: Accuracy test mode by comparing Hall_Conductivity calculation using Kubo-Streda formula and Berry-Phase.
// 3: Open-mp parallel test.
// 4: Free form
#define TEST_MODE 1



// **** Input Settings ****
// Set input parameter and domain. (0?:s-Hamiltonian, 1?:d-Hamiltonian, 2?:impurity, 9?:system and environment)
// -- 00:Exchange,           01:Rashba-SOI,        02:In-Plane HOPPING,    03:Out-of-Plene Hopping,
// -- 04:Polar-angle(theta), 05:Azmual-angle(phi), 06:4s Lorentzian width
// -- 10:Hund coupling,      11:Atomic-SOI,        12:d-orbital T2g level, 13:Crystal field splitting (Eg-T2g)
// -- 14:3d Lorentzian width 15:Tetragonal crystal field splitting
// -- 20:Impurity density,   21:Hybridization, 
// -- 90:Fermi energy(Chemical potential),     91:Filling(particle/unit cell):[0.0:2.0],    92:Layer [1:99]
#define INPUT 90
#define INPUT_MIN -1.5
#define INPUT_MAX 3.5
#define INPUT_DIV 50



// **** Output(physical quantity) ****
// -- 0:DC_CONDUCTIVITY, 1:HALL_CONDUCTIVITY, 2:SPIN_HALL_CONDUCTIVITY
// -- 3:DOS, 4:MAGNETIZATION, 5:TOTAL_ENERGY, 6:AMR_AMPLITUDE
// -- 7:SpinOrbitTorque, 8:SpinAccumulation, 9:Anisotropic_Energy
#define OUTPUT 3
#define PRECISION 14

// **** Special settings for magnetization direction/rotation ****
// For input (04) and (05): Angle rotation axis
#define ROTATION_AXIS 'z'
// For output (9): Fourfold AMR
#define FOURFOLD 1



// **** Hamiltonian ****
// Set system configlation. (0?: TBSQ_Bilayer, 1?: TBSQ_Trilayer, 5?: TBSQ_Periodic, 8?: Dirac fermion 9?: Electron gas model)
// -- 00: FM/NI,  01: NM/NI, 02: FM/NM, 03: FM/TI, 04: 2DTI
// -- 90: 3DEG,   91: 2DEG,  92: 2DEG + QW
// TMSQ = Tight-binding square lattice, QW = Quantum well
#define H_TYPE 90
#define SD_SCATTERING 0
#define LAYER_MAX 99

// Initial parameters
// Conduction electron system
#define INIT_LAYER 1
#define INIT_IPHOPPING 1.0
#define INIT_OPHOPPING 1.0
#define INIT_EXCHANGE 0.5
#define INIT_SOC 0.0
#define INIT_THETA 0
#define INIT_PHI 90
#define INIT_GAMMA 0.001

// Impurity system
#define INIT_NIMP 0.1
#define INIT_IMP_LEVEL 2
#define INIT_IMP_WIDTH 0.001
#define INIT_IMP_EXCHANGE 0.0
#define INIT_IMP_SOI -0.0
#define INIT_IMP_CF 0.5
#define INIT_IMP_TETRA_CF 0.2
//#define INIT_IMP_CHARGE 5
// Bare charge (Screened charge) Fe:26 (5)
#define INIT_HYB_SD 1.0

// Environment
#define INIT_FERMI_ENERGY 1.13
#define INIT_FILLING 0.05



// **** Numerical Integration ****
// Momentum integral
#define K_MIN (-M_PI)
#define K_MAX (M_PI)
#define K_DIV 10
// Energy integral
#define E_MIN -3.5
#define E_MAX 1.5
#define E_DIV 50
#define EF_ITERATION 2


// **** Option ****
// Halfway Mode about Fermi_Energy
// 0: Off
// 1: Addopting the Default Fermi-Level
#define EF_HALFWAY_MODE 0

// Halfway Mode about Energy
// 0: Off
// 1: Taking into the Surface term only
#define SF_HALFWAY_MODE 0

// Select Spin (or Torque) component
// 0:all, 1:x, 2:y, 3:z
#define MU_DIR 0

// Magnetization rate of sublattice
#define SUB_MAGNETIZATION -1

// Modulate current direction
#define NODAL_CONDUCTION 1

//#define PERTUB






// **** Alert for setting errors 
#if (OUTPUT == 6 || OUTPUT == 9)
    #if ROTATION_AXIS != 'z'
        #error "Please set ROTATION_AXIS to z."
    #elif (INPUT == 4) || (INPUT == 5)
        #warning "Meaningless setting: MR vs angle."
    #endif
#endif

#endif
