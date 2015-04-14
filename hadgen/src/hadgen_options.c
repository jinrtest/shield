#define HADGEN_LIB_INTERNAL
#include "hadgen.h"
#include "hadgen_common_blocks.h"
#include <string.h>
#include <stdio.h>

const char* hadgen_get_particle_name(int type) {
   static char _name[20];
   switch(type) {
      case PARTICLE_Neutron:
         strcpy(_name, "Neutron");        break;
      case PARTICLE_AntiNeutron:
         strcpy(_name, "Anti-neutron");   break;
      case PARTICLE_AntiProton:
         strcpy(_name, "Anti-proton");    break;
      case PARTICLE_Proton:
         strcpy(_name, "Proton");         break;
      case PARTICLE_PIminus:
         strcpy(_name, "Pi-");            break;
      case PARTICLE_PIplus:
         strcpy(_name, "Pi+");            break;
      case PARTICLE_PI0:
         strcpy(_name, "Pi0");            break;
      case PARTICLE_Kminus:
         strcpy(_name, "K-");             break;
      case PARTICLE_Kplus:
         strcpy(_name, "K+");             break;
      case PARTICLE_K0:
         strcpy(_name, "K0");             break;
      case PARTICLE_AntiK0:
         strcpy(_name, "Anti-K0");        break;
      case PARTICLE_Gamma:
         strcpy(_name, "Gamma");          break;
      case PARTICLE_Deuterium:
         strcpy(_name, "Deuterium");      break;
      case PARTICLE_Tritium:
         strcpy(_name, "Tritium");        break;
      case PARTICLE_He3:
         strcpy(_name, "He-3");           break;
      case PARTICLE_Alpha:
         strcpy(_name, "Alpha");          break;
      case PARTICLE_NONE:
      default:
         strcpy(_name, "Unknown");
   }
   return _name;
}

int hadgen_get_pdg_code(int particle_type)
{
//   printf("pdg: jpart = %i", particle_type);
   switch(particle_type) {
      case PARTICLE_Neutron:
         return 2112;
      case PARTICLE_AntiNeutron:
         return -2112;
      case PARTICLE_AntiProton:
         return -2212;
      case PARTICLE_Proton:
         return 2212;
      case PARTICLE_PIminus:
         return -211;
      case PARTICLE_PIplus:
         return 211;
      case PARTICLE_PI0:
         return 111;
      case PARTICLE_Kminus:
         return -321;
      case PARTICLE_Kplus:
         return 321;
      case PARTICLE_K0:
         return 130;       // Long-lived only
      case PARTICLE_AntiK0:
         return -130;         // Ich habe keine frage
      case PARTICLE_Gamma:
         return 22;        
      case PARTICLE_Deuterium:
         return hadgen_get_pdg_code_nuclei(2,1);
      case PARTICLE_Tritium:
         return hadgen_get_pdg_code_nuclei(3,1);
      case PARTICLE_He3:
         return hadgen_get_pdg_code_nuclei(3,2);
      case PARTICLE_Alpha:
         return hadgen_get_pdg_code_nuclei(4,2);
      case PARTICLE_NONE:
      default:
         return 0;
   }
}

#define Z_CODE 10000
#define A_CODE 10

int hadgen_get_pdg_code_nuclei(float _A, float _Z) 
{
//    printf("pdg for %3.0f %3.0f", _A, _Z);
   int A = _A, Z = _Z;
   if (A > 250) return -99901;
   if (Z > 125) return -99902;
   if (Z < -125) return -99903;
   if (A <= 0) return -99904;
   // Weight number should be greater than charge number
   if (A < Z) return -99905;
   // check for ambiguities
   if ((A == 1) && (Z == 1)) return hadgen_get_pdg_code(PARTICLE_Proton);
   int result = 1000000000;
   if (Z > 0)
   result += A_CODE * A + Z_CODE * Z;
   else
   result += A_CODE * A - Z_CODE * Z;
   if (Z < 0) result = -result;  //for antinuclei
   return result;
}

void hadgen_set_randomseed(int seed) {
   random.IXFIRS = seed;
}
void hadgen_set_luxcount(int luxcnt) {
   other.LUXCNT = luxcnt;
}
void hadgen_set_incidentparticle(int jpart) {
   // check for validity TODO
   inreac.JPART = jpart;
}
void hadgen_set_nuclid(int nuclid) {
   // check for validity TODO
   inreac.NUCLID = nuclid;
}
void hadgen_set_energy(float energy) { // input energy in GeV
   if (energy > 0) inreac.TINT = energy * 1000.;
   else inreac.TINT = -energy * 1000.;
}
void hadgen_set_system(int lanti) {
   if (lanti != 1) lanti = 0; antlab.lanti = lanti;
}
void hadgen_set_lantil(int lantil) {
   antil.lantil = lantil;
}
void hadgen_set_lstar(int lstar) {
   if (lstar != 1) lstar = 0;
   debug.LSTAR = lstar;
}
void hadgen_set_lcasc(int lcasc) {
   if (lcasc != 1) lcasc = 0;
   debug.LCASC = lcasc;
}
void hadgen_set_statisticsnum(int nstat) {
   if (nstat < 0) nstat = -nstat;
   if (nstat == 0) nstat = 1;
   other.NSTAT = nstat;
}
void hadgen_set_aprojectile(float aproj) {
   // Check for validity will be performed later
   hiproj.APROJ = aproj;
}
void hadgen_set_zprojectile(float zproj) {
   hiproj.ZPROJ = zproj;
}
      
int hadgen_get_randomseed() {
   return random.IXFIRS;
}
int hadgen_get_luxcount() {
   return other.LUXCNT;
}
int hadgen_get_incidentparticle() {
   return inreac.JPART;
}
int hadgen_get_nuclid() {
   return inreac.NUCLID;
}
float hadgen_get_energy() {
   return inreac.TINT * 0.001; // output energy in GeV
}
int hadgen_get_system() {
   return antlab.lanti;
}
int hadgen_get_lantil() {
   return antil.lantil;
}
int hadgen_get_LSTAR() {
   return debug.LSTAR;
}
int hadgen_get_LCASC() {
   return debug.LCASC;
}
int hadgen_get_statisticsnum() {
   return other.NSTAT;
}
float hadgen_get_aprojectile() {
   return hiproj.APROJ;
}
float hadgen_get_zprojectile() {
   return hiproj.ZPROJ;
}

