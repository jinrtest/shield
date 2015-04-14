// @(#)root/montecarlo/hadgen:$Id: THadgen.h
// Author: Alexander Timofeev   07/07/2011

/*************************************************************************
 * Written by Alexander Timofeev, JINR, Dubna                            *
 * On any questions feel free to contact: antimofeew@gmail.com           *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_Hadgen_header
#define ROOT_Hadgen_header

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// THadgen                                                              //
// is a wrapper class for HADGEN generator                              //
// original code for HADGEN (SHIELD) is                                 //
// written by N.M.Sobolevsky                                            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "TObjArray.h"
#include "TClonesArray.h"
#include "TGenerator.h"
#include "TParticle.h"

class THadgen : public TGenerator {
   ClassDef(THadgen, 1)	// THadgen
protected:
   // auto-cleanup mechanism
   static THadgen *fInstance;
   class THadgenClean {
   public:
      THadgenClean() {}
      ~THadgenClean();
   };
   friend class THadgenClean;

   Int_t ParticlesNum;
   bool autoseed;

private:	
   // THadgen class should be never copied
   THadgen(const THadgen& s) : TGenerator(s) {}
   THadgen& operator=(const THadgen&) { return *this; }

public:
   THadgen();		// creates a unique instance of THadgen
   virtual ~THadgen();
	
   static THadgen * Instance();

   Int_t DoSome();	// Just a test function, will be deleted in future release

   // Random generator state
    void SetRandomSeed(Int_t seed);
    void SetRandomSeed();
    Int_t GetRandomSeed();

    void SetAutoSeed(bool flag = true);   // initially true    

protected:
   // the following functions should not be used
    void SetLuxCount(Int_t luxcnt);
    Int_t GetLuxCount();

public:
   // Inreac properties
    void SetIncidentParticle(Int_t jpart);  // these functions may be more smart 
    Int_t GetIncidentParticle();
    void SetNuclid(Int_t nuclid);
    Int_t GetNuclid();
    void SetEnergy(Float_t energy);
    Float_t GetEnergy();
   // Flags
    void SetSystem(Int_t lanti = 0);
    Int_t GetSystem();
    void SetLAntil(Int_t lantil = 0);
    Int_t GetLAntil();
   // Debug options (prInt_ting options)
   // 0 means no output, 1 - force output
    void SetLSTAR(Int_t lstar = 0);
    Int_t GetLSTAR();
    void SetLCASC(Int_t lcasc = 0);
    Int_t GetLCASC();

   // Monte-Carlo statistics parameter
    void SetNumStat(Int_t nstat = 1);
    Int_t GetNumStat();
 
   // A and Z of projectile ion
   // notice that these params are set automaticaly for 
   // Deuteron, Tritium, He3 and alpha particles
    void SetAProj(Float_t aproj);
    Float_t GetAProj();
    void SetZProj(Float_t zproj);
    Float_t GetZProj();

   // Parameter functions derived from TGenerator
//   virtual void SetParameter(const char *name, Double_t value);
//   virtual Double_t GetParameter(const char *name);


   // Generation Functions derived from TGenerator
   void GenerateEvent();
   TObjArray* ImportParticles(Option_t *option="");
   Int_t ImportParticles(TClonesArray* particles, Option_t* option="");
   
   void FileOut(char *filename); // test function
   
};


#endif
