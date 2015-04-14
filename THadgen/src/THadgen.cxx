// @(#)root/montecarlo/hadgen:$Id: THadgen.cxx
// Authors: Alexander Timofeev 07/07/2011

/*************************************************************************
 * Written by Alexander Timofeev, JINR, Dubna                            *
 * On any questions feel free to contact: antimofeew@gmail.com           *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// THadgen                                                                    //
// is a Interface (wrapper) class for Hadgen generator package                //
// Original code of HADGEN (SHIELD) package                                   //
// is written by N.M.Sobolevsky                                               //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "THadgen.h"
#include "hadgen.h"     // link to $SIMPATH/generators/hadgen/hadgen.h
                        // includes all the lib Interface
#include <cstdlib>
#include <time.h>

// class THadgen implementation
THadgen *THadgen::fInstance = 0;


//______________________________________________________________________________
THadgen::THadgenClean::~THadgenClean() {
   // delete the unique copy of THadgen
   // if exists
   if (THadgen::fInstance)
      delete THadgen::fInstance;
}
//______________________________________________________________________________


ClassImp(THadgen)
//______________________________________________________________________________
THadgen::THadgen() : TGenerator("THadgen", "THadgen"), autoseed(true)
{
   // Check if another instance of THadgen already exists
   if (fInstance)
      Fatal("THadgen", "FATAL ERROR: Another instance of THadgen class already exists");

   // init the random generator
   srand(time(0));
   fInstance = this;
   // create static cleaner
   static THadgenClean cleaner;
   // Some initialization of base class needed
   delete fParticles;
   fParticles = new TClonesArray("TParticle", 
      hadgen_get_max_stp() + hadgen_get_max_snu()); 

   // Setting all default
   hadgen_set_defaults();
}

THadgen* THadgen::Instance() {
   return fInstance ? fInstance : new THadgen();
}

//______________________________________________________________________________
THadgen::~THadgen() 
{
   // if deleted, clear the instance
   THadgen::fInstance = 0;
   // useful if needed to create a new instance of THadgen in future
   // actually THadgen object should be automatically cleaned
   // but user can also delete it manually
}

#include <stdio.h>

//______________________________________________________________________________
void THadgen::GenerateEvent() {
   if (autoseed) SetRandomSeed();
/*   if (hadgen_start()) {
      printf("THadgen: initialization error, bad parameters on input.\n");
   } else ImportParticles();*/
   // This is a previous version of hadgen call
     // Actually should work, but use hadgen_start() in case if the generator hangs
   if (hadgen_initialize()) {
      //report error
   } else {
      hadgen_generate();
      hadgen_terminate();
      // copy all data into array of THadgen class 
      ImportParticles();
   }

}

//______________________________________________________________________________
TObjArray *THadgen::ImportParticles(Option_t *) {
   const Int_t MAX_PARTICLES = hadgen_get_max_stp();
   const Int_t MAX_NUCLEI    = hadgen_get_max_snu();

   fParticles->Clear();
   TClonesArray &a = *((TClonesArray*)fParticles);
   Int_t n_spt = 0;
      
   {
      hadgen_iter_stp_reset();   // reset the iterator
      struct HadgenParticle_t buf;
      for(Int_t i=0; i<MAX_PARTICLES; i++) {
         hadgen_iter_stp(&buf);
         if (buf.type != PARTICLE_NONE) {
            Int_t pdg = hadgen_get_pdg_code(buf.type);
//          Debug printout
//            printf("THadgen: imported particle pdg: %i type: %i\n", pdg, buf.type);
            new (a[i]) TParticle(pdg, 0, 0, 0, 0, 0, 
                                 buf.Px, buf.Py, buf.Pz, buf.Energy + buf.Weight,   // TOTAL ENERGY STORED
                                 0, 0, 0, 0);
            n_spt ++;
         } else break;
      }
   }
   Int_t n_snu = 0;
   {
      hadgen_iter_snu_reset();   // reset second iterator for nuclei
      struct HadgenResidualNuclei_t buf;
      for(Int_t i=0; i<MAX_NUCLEI; i++) {
         hadgen_iter_snu(&buf);
         if (buf.A != 0) {
            Int_t pdg = hadgen_get_pdg_code_nuclei(buf.A, buf.Z);
//          Debug printout
//            printf("THadgen: imported nuclei pdg: %i A: %f Z: %f\n", pdg, buf.A, buf.Z);
            new (a[i+n_spt]) TParticle(pdg, 0, 0, 0, 0, 0,
                                       buf.Px, buf.Py, buf.Pz, buf.Energy + 0.940*buf.A, // TOTAL ENERGY
                                       0, 0, 0, 0);
            n_snu ++;
         } else break;
      }
   }
   ParticlesNum=n_spt + n_snu;
   return fParticles;
}

//______________________________________________________________________________
Int_t THadgen::ImportParticles(TClonesArray *particles, Option_t *)
{
   TClonesArray &a = *particles;
   const Int_t MAX_PARTICLES = hadgen_get_max_stp();
   const Int_t MAX_NUCLEI    = hadgen_get_max_snu();

   Int_t n_spt = 0;
   {
      hadgen_iter_stp_reset();   // reset the iterator
      struct HadgenParticle_t buf;
      for(Int_t i=0; i<MAX_PARTICLES; i++) {
         hadgen_iter_stp(&buf);
         if (buf.type != PARTICLE_NONE) {
            Int_t pdg = hadgen_get_pdg_code(buf.type);
            new (a[i]) TParticle(pdg, 0, 0, 0, 0, 0, 
                                 buf.Px, buf.Py, buf.Pz, buf.Energy + buf.Weight,
                                 0, 0, 0, 0);
            n_spt ++;
         } else break;
      }
   }
   Int_t n_snu = 0;
   {
      hadgen_iter_snu_reset();   // reset second iterator for nuclei
      struct HadgenResidualNuclei_t buf;
      for(Int_t i=0; i<MAX_NUCLEI; i++) {
         hadgen_iter_snu(&buf);
         if (buf.A != 0) {
            Int_t pdg = hadgen_get_pdg_code_nuclei(buf.A, buf.Z);
            new (a[i+n_spt]) TParticle(pdg, 0, 0, 0, 0, 0,
                                       buf.Px, buf.Py, buf.Pz, buf.Energy + 0.940*buf.A,
                                       0, 0, 0, 0);
            n_snu ++;
         } else break;
      }
   }

   return n_spt + n_snu;
}

//______________________________________________________________________________
void THadgen::FileOut(char *filename) {
   FILE *f;
   f = fopen(filename, "w");
   fprintf(f, "Number of particles = %i\n", ParticlesNum);
   TClonesArray &a = *((TClonesArray*)fParticles);
   for(Int_t i=0; i<ParticlesNum; i++) {
      TParticle &p = (*(TParticle*)a[i]);
      fprintf(f, "%10i %14e %14e %14e %14e\n", p.GetPdgCode(), 
         p.Px(), p.Py(), p.Pz(), p.Energy());
   }
   fclose(f);   
}

