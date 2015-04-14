#define HADGEN_LIB_INTERNAL

#include <string.h>
#include "hadgen_common_blocks.h"
#include "hadgen.h"
// SOME STUFF USED FOR HADGEN SPECTRA OUTPUT

void Name(int JPART, char buff[]) 
{
   switch(JPART) {
      case 1: // remember all indices start from 1 (FORTRAN)
         strcpy(buff, "Neutron"); break;
      case 2:
         strcpy(buff, "Proton"); break;
      case 3:
         strcpy(buff, "PI-minus"); break;
      case 4:
         strcpy(buff, "PI-plus"); break;
      case 5:
         strcpy(buff, "PI-zero"); break;
      case 6:
         strcpy(buff, "Antineutron"); break;
      case 7:
         strcpy(buff, "Antiproton"); break;
      case 8:
         strcpy(buff, "K-minus"); break;
      case 9:
         strcpy(buff, "K-plus"); break;
      case 10:
         strcpy(buff, "K0"); break;
      case 11:
         strcpy(buff, "AntiK0"); break;
      default:
         strcpy(buff, "UNKNOWN");
   }
}

void hadgen_full_out(struct HadgenOutput_t *result) 
// This dunction fully substitutes the FORTRAN function PRIHAD
{
//   HadgenOutputHad_t result;
   result->JPART = inreac.JPART;
   result->Energy = inreac.TINT;
   if (result->JPART <= 11) {
      result->Antilab = antlab.lanti;
      Name(result->JPART, result->Name);
   } else {
      result->IAPROJ = hiproj.APROJ;
      result->IZPROJ = hiproj.ZPROJ;
      strcpy(result->Name, infel.SYMB[result->IZPROJ-1]);
   }

   strcpy(result->NuclidName, infel.SYMB[inreac.NUCLID-1]);
   result->Znuclid = infel.ZNUC[inreac.NUCLID-1];
   result->Atwei = infel.ATWEI[inreac.NUCLID-1];
   
   // Calculation of cross sections
   float SGEOM, APR, ZPR, ATG, ZTG;
   sigeom(&SGEOM, &APR, &ZPR, &ATG, &ZTG);
   if ((APR <= 1.1) && (ATG <= 1.1)) {
      result->SITO = hhsigm.SITO;
      result->SIEL = hhsigm.SIEL;
      result->SIIN = result->SITO - result->SIEL;
   } else {
      result->INTEL = numint.INTEL;
      result->INTIN = numint.INTIN;
      result->INTOT = result->INTIN + result->INTEL;
   }
   // copying all spectra
   int j;
   for(j=0; j<20; j++) {
      result->DTN[j]  = hadser.DTN[j];
      result->DTP[j]  = hadser.DTP[j];
   }
   for(j=0; j<28; j++) {
      result->SN05[j] = hadron.SN05[j];
      result->SN30[j] = hadron.SN30[j];
      result->SN90[j] = hadron.SN90[j];
      result->SN18[j] = hadron.SN18[j];
      result->SP05[j] = hadron.SP05[j];
      result->SP30[j] = hadron.SP30[j];
      result->SP90[j] = hadron.SP90[j];
      result->SP18[j] = hadron.SP18[j];
      result->SG18[j] = hadron.SG18[j];
      result->PM05[j] = hadron.PM05[j];
      result->PM30[j] = hadron.PM30[j];
      result->PM90[j] = hadron.PM90[j];
      result->PM18[j] = hadron.PM18[j];
      result->PP05[j] = hadron.PP05[j];
      result->PP30[j] = hadron.PP30[j];
      result->PP90[j] = hadron.PP90[j];
      result->PP18[j] = hadron.PP18[j];
      result->P018[j] = hadron.P018[j];
      result->SAN[j]  = hadron.SAN[j];
      result->SAP[j]  = hadron.SAP[j];
      result->SKM[j]  = hadron.SKM[j];
      result->SKP[j]  = hadron.SKP[j];
      result->SK0[j]  = hadron.SK0[j];
      result->SAK0[j] = hadron.SAK0[j];
      result->SLNT[j] = hadron.SLNT[j];
      result->SLPR[j] = hadron.SLPR[j];
      result->SLDT[j] = hadron.SLDT[j];
      result->SLTR[j] = hadron.SLTR[j];
      result->SLH3[j] = hadron.SLH3[j];
      result->SLAL[j] = hadron.SLAL[j];
   }

   result->IFISS = hadron.IFISS;
   result->IDES  = hadron.IDES;

   // Excitation and kinetic energy spectra
   for(j=0; j<28; j++) {
      result->USTR[j]  = hadron.USTR[j];
      result->USTR2[j] = hadron.USTR2[j];
      result->USTR5[j] = hadron.USTR5[j];
      result->TSTR[j]  = hadron.TSTR[j];
      result->TSTRH[j] = hadron.TSTRH[j];
   }
   // Residual nuclei (A,Z)-distribution
   // normally the end-user should scan 
   // all arrays in the range IMIN0 <= I < IMAX0
   result->IMIN0 = 12;
   result->IMAX0 = 12; 
   int IMIN=0, IMAX=24;
   for(j=0; j<100; j++) {
      for(; ((IMIN<25)   && (hadron.IZA[j][IMIN] == 0)); IMIN++);
      for(; ((IMAX>=0)   && (hadron.IZA[j][IMAX] == 0)); IMAX--);
      if (result->IMIN0 > IMIN) result->IMIN0 = IMIN;
      if (result->IMAX0 < IMAX) result->IMAX0 = IMAX;
   }
   // copying all data
   int i=0;
   for(j=0; j<25; j++) result->IAA[j] = resn01.IAA[j];
   for(j=0; j<100; j++) {
      for(i=0; i<100; i++) result->IZA[j][i] = hadron.IZA[j][i];
   }
// PRICON output not implemented YET

}

