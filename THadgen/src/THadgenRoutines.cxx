// @(#)root/montecarlo/hadgen:$Id: THadgen.cxx
// Authors: Alexander Timofeev 07/07/2011

/*************************************************************************
 * Written by Alexander Timofeev, JINR,  Dubna                           *
 * On any questions feel free to contact: antimofeew@gmail.com           *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// THadgen                                                                    //
// is an interface (wrapper) class for Hadgen generator package               //
// Original code for HADGEN (SHIELD) package                                  //
// is written by N.M.Sobolevsky                                               //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "THadgen.h"

#include "hadgen.h" 	// link to $SIMPATH/generators/hadgen/hadgen.h
			            // includes all lib interface
#include <cstdlib>
#include <time.h>

// setters in THadgen class
//______________________________________________________________________________
void THadgen::SetRandomSeed(Int_t seed) {
   if (seed < 0) seed = -seed;   //should be positive!
   hadgen_set_randomseed(seed);
}
//______________________________________________________________________________
void THadgen::SetRandomSeed() {
   SetRandomSeed(clock()+time(0));
}
//______________________________________________________________________________
void THadgen::SetAutoSeed(bool flag) {
   autoseed = flag;
}
//______________________________________________________________________________
void THadgen::SetLuxCount(Int_t luxcnt) {
   hadgen_set_luxcount(luxcnt);
}
//______________________________________________________________________________
void THadgen::SetIncidentParticle(Int_t jpart) {
   hadgen_set_incidentparticle(jpart);
}
//______________________________________________________________________________
void THadgen::SetNuclid(Int_t nuclid) {
   hadgen_set_nuclid(nuclid);
}
//______________________________________________________________________________
void THadgen::SetEnergy(Float_t energy) {
   hadgen_set_energy(energy);
}
//______________________________________________________________________________
void THadgen::SetSystem(Int_t lanti) {
   hadgen_set_system(lanti);
}
//______________________________________________________________________________
void THadgenSetLAntil(Int_t lantil) {
   hadgen_set_lantil(lantil);
}
//______________________________________________________________________________
void THadgen::SetLSTAR(Int_t lstar) {
   hadgen_set_lstar(lstar);
}
//______________________________________________________________________________
void THadgen::SetLCASC(Int_t lcasc) {
   hadgen_set_lcasc(lcasc);
}
//______________________________________________________________________________
void THadgen::SetNumStat(Int_t nstat) {
   hadgen_set_statisticsnum(nstat);
}
//______________________________________________________________________________
void THadgen::SetAProj(Float_t aproj) {
   hadgen_set_aprojectile(aproj);
}
//______________________________________________________________________________
void THadgen::SetZProj(Float_t zproj) {
   hadgen_set_zprojectile(zproj);
}

// getters in THadgen
//______________________________________________________________________________
Int_t THadgen::GetRandomSeed() {
   return hadgen_get_randomseed();
}
//______________________________________________________________________________
Int_t THadgen::GetLuxCount() {
   return hadgen_get_luxcount();
}
//______________________________________________________________________________
Int_t THadgen::GetIncidentParticle() {
   return hadgen_get_incidentparticle();
}
//______________________________________________________________________________
Int_t THadgen::GetNuclid() {
   return hadgen_get_nuclid();
}
//______________________________________________________________________________
Float_t THadgen::GetEnergy() {
   return hadgen_get_energy();
}
//______________________________________________________________________________
Int_t THadgen::GetSystem() {
   return hadgen_get_system();
}
//______________________________________________________________________________
Int_t THadgen::GetLAntil() {
   return hadgen_get_lantil();
}
//______________________________________________________________________________
Int_t THadgen::GetLSTAR() {
   return hadgen_get_lstar();
}
//______________________________________________________________________________
Int_t THadgen::GetLCASC() {
   return hadgen_get_lcasc();
}
//______________________________________________________________________________
Int_t THadgen::GetNumStat() {
   return hadgen_get_statisticsnum();
}
//______________________________________________________________________________
Float_t THadgen::GetAProj() {
   return hadgen_get_aprojectile();
}
//______________________________________________________________________________
Float_t THadgen::GetZProj() {
   return hadgen_get_zprojectile();
}
