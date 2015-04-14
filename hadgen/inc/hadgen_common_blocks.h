// HADGEN library reimplementation header
// needed to connect HADGEN to CERN FAIRROOT
// Written by Alexander Timofeev, 08/07/2011
// Copyright (c) 2011 Alexander Timofeev, Dubna

#ifndef HADGEN_COMMON_BLOCKS_HEADER
#define HADGEN_COMMON_BLOCKS_HEADER

// Definition of HADGEN common blocks
#ifndef WIN32
// Linux implementation
#define random       random_
#define sechar       sechar_
#define inreac       inreac_
#define specagt      specagt_
#define numint       numint_
#define hiproj       hiproj_
#define antlab       antlab_
#define antil        antil_
#define debug        debug_
#define islerr       islerr_
#define other        other_
#define hadron       hadron_
#define resn01       resn01_
#define infel        infel_
#define hhsigm       hhsigm_
#define hadser       hadser_
#define abunuc       abunuc_
#define atabdata     atabdata_
#else
// WIN32
#define random       RANDOM
#define sechar       SECHAR
#define inreac       INREAC
#define specagt      SPECAGT
#define numint       NUMINT
#define hiproj       HIPROJ
#define antlab       ANTLAB
#define antil        ANTIL
#define debug        DEBUG
#define islerr       ISLERR
#define other        OTHER
#define hadron       HADRON
#define resn01       RESN01
#define infel        INFEL
#define hhsigm       HHSIGM
#define hadser       HADSER
#define abunuc       ABUNUC
#define atabdata     ATABDATA
#endif

// EVERY THING IS IN SINGLE PRECISION
struct Random_t {
   int IX;
   int IXINIT;
   int NUMTRE;
   int IXFIRS;
   int LASTIX;
} extern random;

struct Sechar_t {
   float SPT[5000][6];
   float SNU[101][10];
   int LS6;
   int LS100;
   int LS10;
   int LS11;
} extern sechar;

struct Inreac_t {
   float COST;
   float SINF;
   float COSF;
   float TINT;
   float WINT;
   int JPART;
   int NUCLID;
   int KSTATE;
} extern inreac;

struct Specagt_t {
// TO DO!!
// Fortran definition is SPECATG(126,30,180)
// Should be checked
   float SPECAGT[180][30][126];
} extern specagt;

struct Numint_t {
   int INTIN;
   int INTEL;
} extern numint;

struct Hiproj_t {
   float APROJ;
   float ZPROJ;
} extern hiproj;

struct Antlab_t {
   int lanti;
} extern antlab;

struct Antil_t {
   int lantil;
} extern antil;

struct Debug_t {
   int LSTAR;
   int LCASC;
} extern debug;

struct Islerr_t {
   int ISLERR;
} extern islerr;

struct Other_t {
   int NSTAT, LUXCNT;
} extern other;

struct Hadron_t {
   float SN05[28], SN30[28], SN90[28], SN18[28],
         SP05[28], SP30[28], SP90[28], SP18[28], SG18[28],
         PM05[28], PM30[28], PM90[28], PM18[28],
         PP05[28], PP30[28], PP90[28], PP18[28], P018[28],
         SLNT[28], SLPR[28], SLDT[28], SLTR[28], SLH3[28], SLAL[28],
         SAN[28],  SAP[28],  SKM[28],  SKP[28],  SK0[28],  SAK0[28],
         BUSTR[28], BUSTR2[28], BUSTR5[28], BTSTR[28], BTSTRH[28],
         USTR[28], USTR2[28], USTR5[28], TSTR[28], TSTRH[28];
   int IFISS, IDES,
       IZAB[100][25], IOUTB, IINB,
       IZA[100][25], IOUT, IIN; 
} extern hadron;

struct Hadser_t {
   float TMAXN;
   float HN;
   float TMAXP;
   float HP;
   float DTN[20];
   float DTP[20];   
} extern hadser;

struct Resn01_t {
   float ASTAB[100];
   int   IAA[25];
} extern resn01;

struct Infel_t {
   float ZNUC[110];
   float ATWEI[110];
   float DENS[110];
   float RIEV[110];
   char  SYMB[110][2];	// TODO!!! A part to CHECK
} extern infel;

struct HHSigm_t {
   float SITO, SIEL;
} extern hhsigm;

struct Abunuc_t {
   float AIST[500];
   float PERC[500];
   float ABUNOR[500];
   int NUCADR[110][2];
   int IZMAX;
   int L110;
} extern abunuc;

struct Atabdata_t {
   int IRES[600];
   int ITYPE[600];
   double BR[600];
   int IMODE[600][5];
} extern atabdata;


// internal C-functions added into HADGEN library
#ifndef HADGEN_LIB_INTERNAL
#define HADGEN_API extern
#else
#define HADGEN_API
#endif
HADGEN_API void _hadgen_initis();
HADGEN_API void *hadgen_thread(void *);
HADGEN_API void _hadgen_setdky();

// Definition of FORTRAN functions
#ifndef WIN32
 #define fortran_call
 #define clehad 		         clehad_
 #define cleion 	         	cleion_
 #define initam 	         	initam_
 #define initis	         	initis
 #define inlevr		      	inlevr_
 #define init_lux	         	initlux2_
 #define makegrid	         	makegrid_
 #define rluxat		      	rluxat_
 #define react   	         	react_
 #define sigeom		      	sigeom_
 #define hadgen_standalone 	hadgen_
 #define hadgen_short         hadgenshort_
 #define openfiles            openfiles_
 #define closefiles           closefiles_
#else
 #define fortran_call 	   	_stdcall
 #define clehad               CLEHAD
 #define cleion               CLEION
 #define initam               INITAM 
 #define initis               INITIS 
 #define inlevr               INLEVR 
 #define init_lux	         	INITLUX2
 #define makegrid		         MAKEGRID
 #define rluxat	            RLUXAT 
 #define react                REACT
 #define sigeom               SIGEOM
 #define hadgen_standalone	   HADGEN
 #define hadgen_short         HADGENSHORT
 #define openfiles            OPENFILES
 #define closefiles           CLOSEFILES
#endif


// MAIN HADGEN function
// May be used to run HADGEN as a standalone application
// Refer to original source code for input and output files
extern void fortran_call hadgen_standalone();

extern void fortran_call hadgen_short();

extern void fortran_call clehad();
extern void fortran_call cleion();
extern void fortran_call initam();
extern void fortran_call initis();
extern void fortran_call inlevr();
extern void fortran_call init_lux(int*, int*);
extern void fortran_call makegrid();
extern void fortran_call rluxat(int*, int*, int*, int*);
extern void fortran_call react(int *);
extern void fortran_call sigeom(float*, float*, float*, float*, float*);

extern void fortran_call openfiles();
extern void fortran_call closefiles();

// NOTICE
// All spectra contain two parts
// points 1-20 (fortran index) contain energy distribution
// points 21-28 (fortran index) contain other data

// the struct for output
// all data should be COPIED from the fortran code, NOT LINKED via pointers

#endif
