// File:    shield_common_blocks.h
// Authors: Alexander Timofeev 25/12/2012
//          Dmitry    Sosnov      02/2015

#ifndef SHIELD_COMMON_BLOCKS_H
#define	SHIELD_COMMON_BLOCKS_H

#include "shield.h"
#include "hadgen_common_blocks.h"
// this dependency should be removed in future


#ifndef oln
#define oln oln_
struct Oln_t {
    float OLN;
} extern oln;
#endif

#ifndef t0jpr0
#define t0jpr0 t0jpr0_
struct T0jpr0_t {
    float TMAX0;
    int JPART0;
} t0jpr0;
#endif

#ifndef isoset
#define isoset isoset_
struct Isoset_t {
    float WISOSET[2000];
    float FERTAZ[50][2000];
    int ISOSET[2000][2];
    char TDEC[2000][10];
    int NIS200;
    int NISTRU;
} extern isoset;
#endif

#ifndef inre00
#define inre00 inre00_
extern struct Inreac_t inre00;
#endif

#ifndef hipro0
#define hipro0 hipro0_
extern struct Hiproj_t hipro0;
#endif


// ************** Shield result output ***************
// COMMON-blocks

#ifndef tree
#define tree tree_
struct Tree_t {
    float TREE[400000][15];
    int ITREE[400000][7];
    int LTR13;
    int LTR7;
    int LTR300;
} extern tree;
#endif

#ifndef sourn
#define sourn sourn_
struct Sourn_t {
    float SOURN[100000][13];
    int ISOURN[100000][7];
    int LSN13;
    int LSN7;
    int LSN100;
} extern sourn;
#endif

#ifndef sourg
#define sourg sourg_
struct Sourg_t {
    float SOURG[20000][13];
    int ISOURG[20000][7];
    int LSG13;
    int LSG7;
    int LSG100;
} extern sourg;
#endif

#ifndef semin 
#define semin semin_
struct Semin_t {
    float SEMIN[100][13];
    int ISEMIN[100][7];
    int L1313;
    int L137;
    int L13100;
} extern semin;
#endif

#ifndef sepls
#define sepls sepls_
struct Sepls_t {
    float SEPLS[100][13];
    int ISEPLS[100][7];
    int L1413;
    int L147;
    int L14100;
} extern sepls;
#endif

#ifndef snuel
#define snuel snuel_
struct Snuel_t {
    float SNUEL[100][13];
    int ISNUEL[100][7];
    int L1713;
    int L177;
    int L17100;
} extern snuel;
#endif

#ifndef sanue
#define sanue sanue_
struct Sanue_t {
    float SANUE[100][13];
    int ISANUE[100][7];
    int L1813;
    int L187;
    int L18100;
} extern sanue;
#endif

#ifndef snumu
#define snumu snumu_
struct Snumu_t {
    float SNUMU[100][13];
    int ISNUMU[100][7];
    int L1913;
    int L197;
    int L19100;
} extern snumu;
#endif


#ifndef sanum
#define sanum sanum_
struct Sanum_t {
    float SANUM[100][13];
    int ISANUM[100][7];
    int L2013;
    int L207;
    int L20100;
} extern sanum;
#endif

#ifndef trevs
#define trevs trevs_
struct Trevs_t {
    float TREVS[400][15];
    int ITREVS[400][7];
    int LVS13;
    int LVS7;
    int LVS100;
} extern trevs;
#endif

// Initial position common-block
#ifndef wrpr
#define wrpr wrprin_
struct Wrpr_t {
    float WP[15];
    int IWP[7];
} extern wrpr;
#endif

// Media COMMON-block
#ifndef macrin
#define macrin macrin_
struct Macrin_t {
    float ELMIX[48][24][8];
    float RHOMED[48];
    int NELEMD[48][2];
    int NUMMED;
} extern macrin;
#endif

// Media data
#include "shield_media.h"
extern int MediaNum;

// Geometry block
#ifndef gdata0
#define gdata0 gdata0_
struct Gdata0_t {
    int NUMBOD;
    int NUMZON;
} extern gdata0;
#endif
#ifndef gdata2
#define gdata2 gdata2_
struct Gdata2_t {
    int BODYDA[1000];
    double BODYDB[37000];
} extern gdata2;
#endif
#ifndef gdata3
#define gdata3 gdata3_
struct Gdata3_t {
    int ZONEDA[1000];
    int ZONEDB[5000];
} extern gdata3;
#endif

// internal C functions
void shield_set_iso();
void shield_add_iso();
void shield_init_medium();
void shield_init_geometry();
void shield_clean_medium(void);

// internal FORTRAN functions
#ifndef WIN32
#define fortran_call
#define dimar  dimar_
#define calmac calmac_
#define roptpr roptpr_
#define bdprin bdprin_
#define calmhi calmhi_
#define ropthi ropthi_
#define inmacr inmacr_
#define lnmacr lnmacr_
#define cletot cletot_
#define clerln clerln_
#define geoini geoini_
#define cutoff cutoff_
#define gentre gentre_
#define outtot outtot_
#define neutlo neutlo_
#define dectre dectre_
//     #define gcurzl gcurzl_      // GEMCA subsitution
//     #define gnextz gnextz_      // GEMCA subsitution
//     #define gzmed  gzmed_      // GEMCA subsitution
#define opensh opensh_
#define clossh clossh_
#define init_luxs initlux3_
#else
#define fortran_call __stdcall
#endif

extern void fortran_call dimar();
extern void fortran_call calmac();
extern void fortran_call roptpr();
extern void fortran_call bdprin();
extern void fortran_call calmhi(int*);
extern void fortran_call ropthi(int*);
extern void fortran_call lnmacr();
extern void fortran_call cletot();
extern void fortran_call clerln();
extern void fortran_call geoini();      // GEMCA initialization, should be changed
extern void fortran_call cutoff();
extern void fortran_call inmacr();
extern void fortran_call gentre(int*);
extern void fortran_call outtot(int*, float*);
extern void fortran_call neutlo();
extern void fortran_call dectre();
extern void fortran_call init_luxs(int *, int *);

extern void fortran_call opensh();
extern void fortran_call clossh();

#ifndef SHIELD_TREE_DECL
#define SHIELD_TREE_DECL
#define shield_tree_reset       treeres_
#define shield_tree_new_branch  treedec_
void shield_tree_new_branch(int *kstate, int *j_start, int *j_end, int *n_part_cas, int *n_part_tot);
void shield_tree_reset();
#endif

// Geometry module
// Instead of GEMCA
// TODO test if double* or float* is used
void gcurzl(int *, float*); // current zone length
void gnextz(int *);         // next zone number
void gzmed(int *);

#endif	/* SHIELD_COMMON_BLOCKS_H */

