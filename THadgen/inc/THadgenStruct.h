// @(#)root/montecarlo/hadgen:$Id: THadgenStruct.h
// Author: Alexander Timofeev   06/07/2011

/*************************************************************************
 * Written by Alexander Timofeev, JINR, Dubna                            *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_HadgenStruct_header
#define ROOT_HadgenStruct_header

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// THadgen                                                              //
// NOTICE: This header is deprecated                                    //
//////////////////////////////////////////////////////////////////////////

struct Random_t {
   int IX;
   int IXINIT;
   int NUMTRE;
   int IXFIRS;
   int LASTIX;
};

struct Sechar_t {
   float SPT[5000][6];
   float SNU[101][10];
   int LS6;
   int LS100;
   int LS10;
   int LS11;
};

struct Inreac_t {
   float COST;
   float SINF;
   float TINT;
   float WINT;
   int JPART;
   int NUCLID;
   int KSTATE;
};

struct Specagt_t {
// TO DO!!
// Fortran definition is SPECATG(126,30,180)
// Should be checked
   float SPECAGT[180][30][126];
};

struct Numint_t {
   int INTIN;
   int INTEL;
};

struct Hiproj_t {
   float APROJ;
   float ZPROJ;
};

struct Antlab_t {
   int lanti;
};

struct Antil_t {
   int lantil;
};

struct Debug_t {
   int LSTAR;
   int LCASC;
};

struct Islerr_t {
   int ISLERR;
};

struct Other_t {
   int NSTAT, LUXCNT;
};

#endif
