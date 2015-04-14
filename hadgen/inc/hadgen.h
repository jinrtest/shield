// Written by Alexander Timofeev on 12-07-2011, Dubna

#ifdef HADGEN_LIB_INTERNAL
#ifdef __cplusplus
   #error Do not compile HADGEN lib with C++ compiler
#endif
#endif

// INCLUDE FILE FOR USAGE OF HADGEN GENERATOR
#ifndef HADGEN_HEADER
#define HADGEN_HEADER

// Particle types definition according to HADGEN definitions
#define PARTICLE_NONE            0

#define PARTICLE_Neutron         1
#define PARTICLE_Proton          2

#define PARTICLE_PIminus         3
#define PARTICLE_PIplus          4
#define PARTICLE_PI0             5

#define PARTICLE_AntiNeutron     6
#define PARTICLE_AntiProton      7

#define PARTICLE_Kminus          8
#define PARTICLE_Kplus           9
#define PARTICLE_K0              10
#define PARTICLE_AntiK0          11

#define PARTICLE_Gamma           12
#define PARTICLE_Deuterium       21
#define PARTICLE_Tritium         22
#define PARTICLE_He3             23
#define PARTICLE_Alpha           24
#define PARTICLE_Heavy           25

// HadGen particle output structure
struct HadgenParticle_t {
   int type;            // SEE ABOVE
   float Px, Py, Pz;    // In [GeV/c]
   float Energy;        // In [GeV]
   float Charge;        // In charge units
   float Weight;        // In [GeV/c^2]
   float Lifetime;      // In [mm/c]
};

struct HadgenResidualNuclei_t {
   float A;
   float Z;
   float ExcEnergy;     // Excitation energy [GeV]
   float kinter;
   float Px, Py, Pz;    // [GeV/c]
   float Energy;        // Kynetic energy [GeV]
};

// HadGen statistics output structure
struct HadgenOutput_t {
// generic information
   // Incident particle info
   int JPART;		// particle type

   char Name[12]; 	// incident particle (or ion) name NAME(JAPRT) or SYMB(IZPROJ)
   // if JPART <= 11
   int Antilab;		// 1 if Particle system, 0 if Lab. system
   // otherwise (for ions only)
   int IAPROJ, IZPROJ;
   float Energy;	// TINT, particle energy

   // Nucleus target
   char NuclidName[4];	// Nuclid name
   int Znuclid;		// nuclid charge
   float Atwei;		// don't know
   
// cross-sections
   // hh-interaction
   float SIIN, SITO, SIEL;	// Inelastic, total and elastic cross-sections
   // hA- and AA-interactions
   float INTOT, INTIN, INTEL;	// Total, inelastic and elastic cross-sections

// particles spectra   
   // Nucleon and gamma spectra
   float DTN[20], SN05[28], SN30[28], SN90[28], SN18[28], // Neutrons
                  SP05[28], SP30[28], SP90[28], SP18[28], // Protons
                  SG18[28];                               // Gamma
   // Pion spectra
   float DTP[20], PM05[28], PM30[28], PM90[28], PM18[28], // PI-
                  PP05[28], PP30[28], PP90[28], PP18[28], // PI+
                  P018[28];
   // Antineutrons, antiprotons, kaons, antikaons spectra   
   float SAN[28], SAP[28],  // Antineutrons, antiprotons
                  SKM[28], SKP[28],  // K-, K+
                  SK0[28], SAK0[28]; // K0, anti-K0
   // Low energy particle spectra (neutrons, protons, deutrons, tritons, HE3, alpha)
   float SLNT[28], SLPR[28], SLDT[28], SLTR[28], SLH3[28], SLAL[28];

// After deexcitation:
   // Number of fissions and number of full desintegrations   
   int IFISS, IDES;

// RESNPR(USTR,USTR2,USTR5,TSTR,TSTRH,IZA,IOUT,IIN)
   // All spectra are presented in different scales
   // should be reduced to one spectrum
   // TODO!!!
   // Excitation spectrum
   float USTR[28], USTR2[28], USTR5[28]; //0-10, 0-200, 0-500 MeV
   // Kinetic energy spectrum
   float TSTR[28], TSTRH[28];            //0-20, 0-100  MeV
   
// Residual nuclei (A,Z)-ditribution
   // Num of nuclei inside and outside the distribution
   int IIN, IOUT;
   int IMIN0, IMAX0;
   int IAA[25];
   int IZA[100][25];
   
// PRICON subroutine data output
   // NOT IMPLEMENTED YET
};

#ifndef HADGEN_LIB_INTERNAL
    #ifdef HADGEN_DOUBLE_UNDERSCORE
      #define hadgen_get_block               hadgen_get_block__
      #define hadgen_check_params            hadgen_check_params__   
      #define hadgen_initialize              hadgen_initialize__
      #define hadgen_generate                hadgen_generate__
      #define hadgen_terminate               hadgen_terminate__
      #define hadgen_full_out                hadgen_full_out__
      #define hadgen_iter_stp                hadgen_iter_stp__
      #define hadgen_iter_stp_reset          hadgen_iter_stp_reset__
      #define hadgen_get_max_stp             hadgen_get_max_stp__
      #define hadgen_iter_snu                hadgen_iter_snu__
      #define hadgen_iter_snu_reset          hadgen_iter_snu_reset__
      #define hadgen_get_max_snu             hadgen_get_max_snu__
      #define hadgen_get_particle_name       hadgen_get_particle_name__
      #define hadgen_get_pdg_code            hadgen_get_pdg_code__
      #define hadgen_get_pdg_code_nuclei     hadgen_get_pdg_code_nuclei__
      #define hadgen_set_defaults            hadgen_set_defaults__
      #define hadgen_set_random_seed         hadgen_set_random_seed__
      #define hadgen_set_luxcount            hadgen_set_luxcount__
      #define hadgen_set_incidentparticle    hadgen_set_incidentparticle__
      #define hadgen_set_nuclid              hadgen_set_nunclid__
      #define hadgen_set_energy              hadgen_set_energy__
      #define hadgen_set_system              hadgen_set_system__
      #define hadgen_set_lantil              hadgen_set_lantil__
      #define hadgen_set_lstar               hadgen_set_lstar__
      #define hadgen_set_lcasc               hadgen_set_lcasc__
      #define hadgen_set_statisticsnum       hadgen_set_statisticsnum__
      #define hadgen_set_aprojectile         hadgen_set_aprojectile__
      #define hadgen_set_zprojectile         hadgen_set_zprojectile__
      #define hadgen_get_random_seed         hadgen_get_random_seed__
      #define hadgen_get_luxcount            hadgen_get_luxcount__
      #define hadgen_get_incidentparticle    hadgen_get_incidentparticle__
      #define hadgen_get_nuclid              hadgen_get_nunclid__
      #define hadgen_get_energy              hadgen_get_energy__
      #define hadgen_get_system              hadgen_get_system__
      #define hadgen_get_lantil              hadgen_get_lantil__
      #define hadgen_get_lstar               hadgen_get_lstar__
      #define hadgen_get_lcasc               hadgen_get_lcasc__
      #define hadgen_get_statisticsnum       hadgen_get_statisticsnum__
      #define hadgen_get_aprojectile         hadgen_get_aprojectile__
      #define hadgen_get_zprojectile         hadgen_get_zprojectile__
      #define hadgen_print_setup             hadgen_print_setup__
    #elif HADGEN_SINGLE_UNDERSCORE
      #define hadgen_get_block               hadgen_get_block_
      #define hadgen_check_params            hadgen_check_params_   
      #define hadgen_initialize              hadgen_initialize_
      #define hadgen_generate                hadgen_generate_
      #define hadgen_terminate               hadgen_terminate_
      #define hadgen_full_out                hadgen_full_out_
      #define hadgen_iter_stp                hadgen_iter_stp_
      #define hadgen_iter_stp_reset          hadgen_iter_stp_reset_
      #define hadgen_get_max_stp             hadgen_get_max_stp_
      #define hadgen_iter_snu                hadgen_iter_snu_
      #define hadgen_iter_snu_reset          hadgen_iter_snu_reset_
      #define hadgen_get_max_snu             hadgen_get_max_snu_
      #define hadgen_get_particle_name       hadgen_get_particle_name_
      #define hadgen_get_pdg_code            hadgen_get_pdg_code_
      #define hadgen_get_pdg_code_nuclei     hadgen_get_pdg_code_nuclei_
      #define hadgen_set_defaults            hadgen_set_defaults_
      #define hadgen_set_random_seed         hadgen_set_random_seed_
      #define hadgen_set_luxcount            hadgen_set_luxcount_
      #define hadgen_set_incidentparticle    hadgen_set_incidentparticle_
      #define hadgen_set_nuclid              hadgen_set_nunclid_
      #define hadgen_set_energy              hadgen_set_energy_
      #define hadgen_set_system              hadgen_set_system_
      #define hadgen_set_lantil              hadgen_set_lantil_
      #define hadgen_set_lstar               hadgen_set_lstar_
      #define hadgen_set_lcasc               hadgen_set_lcasc_
      #define hadgen_set_statisticsnum       hadgen_set_statisticsnum_
      #define hadgen_set_aprojectile         hadgen_set_aprojectile_
      #define hadgen_set_zprojectile         hadgen_set_zprojectile_
      #define hadgen_get_random_seed         hadgen_get_random_seed_
      #define hadgen_get_luxcount            hadgen_get_luxcount_
      #define hadgen_get_incidentparticle    hadgen_get_incidentparticle_
      #define hadgen_get_nuclid              hadgen_get_nunclid_
      #define hadgen_get_energy              hadgen_get_energy_
      #define hadgen_get_system              hadgen_get_system_
      #define hadgen_get_lantil              hadgen_get_lantil_
      #define hadgen_get_lstar               hadgen_get_lstar_
      #define hadgen_get_lcasc               hadgen_get_lcasc_
      #define hadgen_get_statisticsnum       hadgen_get_statisticsnum_
      #define hadgen_get_aprojectile         hadgen_get_aprojectile_
      #define hadgen_get_zprojectile         hadgen_get_zprojectile_
      #define hadgen_print_setup             hadgen_print_setup_
   #else 
   #endif
   #ifdef __cplusplus
      #define HADGEN_API   extern "C"
   #else
      #define HADGEN_API   extern
   #endif
#else
   // HADGEN INTERNAL
   #define HADGEN_API
#endif

#ifndef HADGEN_API
#define HADGEN_API
#endif

HADGEN_API   int hadgen_start();

HADGEN_API   void *hadgen_get_block(char *block_name); // should not be used actually
HADGEN_API   int hadgen_check_params();
HADGEN_API   int hadgen_initialize(); // return 0 if success, otherwise 1
HADGEN_API   void hadgen_generate();
HADGEN_API   void hadgen_terminate();

HADGEN_API   void hadgen_full_out(struct HadgenOutput_t *);
HADGEN_API   void hadgen_iter_stp(struct HadgenParticle_t *);
HADGEN_API   void hadgen_iter_stp_reset();
HADGEN_API   int hadgen_get_max_stp();

HADGEN_API   void hadgen_iter_snu(struct HadgenResidualNuclei_t *);
HADGEN_API   void hadgen_iter_snu_reset();
HADGEN_API   int hadgen_get_max_snu();

HADGEN_API   const char* hadgen_get_particle_name(int type);
HADGEN_API   int hadgen_get_pdg_code(int particle_type);
HADGEN_API   int hadgen_get_pdg_code_nuclei(float A, float Z);

   // all setters and getters
HADGEN_API   void hadgen_set_defaults();

HADGEN_API   void hadgen_set_randomseed(int seed);
HADGEN_API   void hadgen_set_luxcount(int luxcnt);
HADGEN_API   void hadgen_set_timeout(int sec);
HADGEN_API   void hadgen_set_incidentparticle(int jpart);
HADGEN_API   void hadgen_set_nuclid(int nuclid);
HADGEN_API   void hadgen_set_energy(float energy);
HADGEN_API   void hadgen_set_system(int lanti);
HADGEN_API   void hadgen_set_lantil(int lantil);
HADGEN_API   void hadgen_set_lstar(int lstar); // used only      
HADGEN_API   void hadgen_set_lcasc(int lcasc); // for debug

HADGEN_API   void hadgen_set_statisticsnum(int nstat); 
   // should set 1 for one event, more for statistics
   // initially set to one
   
HADGEN_API   void hadgen_set_aprojectile(float aproj);
HADGEN_API   void hadgen_set_zprojectile(float zproj);
   
HADGEN_API   int hadgen_get_randomseed();
HADGEN_API   int hadgen_get_luxcount();
HADGEN_API   int hadgen_get_timeout();
HADGEN_API   int hadgen_get_incidentparticle();
HADGEN_API   int hadgen_get_nuclid();
HADGEN_API   float hadgen_get_energy();
HADGEN_API   int hadgen_get_system();
HADGEN_API   int hadgen_get_lantil();
HADGEN_API   int hadgen_get_lstar(); // used only      
HADGEN_API   int hadgen_get_lcasc(); // for debug
HADGEN_API   int hadgen_get_statisticsnum();
HADGEN_API   float hadgen_get_aprojectile();
HADGEN_API   float hadgen_get_zprojectile();  

HADGEN_API   void hadgen_print_setup();

#endif // HADGEN_HEADER
