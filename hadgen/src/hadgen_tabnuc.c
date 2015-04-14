#define HADGEN_LIB_INTERNAL
#include "hadgen_common_blocks.h"
#include <stdio.h>

static int isotope_count = 0;
static int Z_count = 0;
static int IAD1 = 0, IAD2 = 0;
static int NISOT = 0;

void _hdg_add_isotope(float Z, float A, float Perc)
{
   abunuc.AIST[isotope_count] = A;
   abunuc.PERC[isotope_count] = Perc;
   NISOT ++;
   isotope_count++;
}

void _hdg_done() {
   IAD1 = IAD2 + 1;
   IAD2 = IAD1 + NISOT - 1;
   NISOT = 0;
   abunuc.NUCADR[Z_count][0] = IAD1;
   abunuc.NUCADR[Z_count][1] = IAD2;
   Z_count++;
   abunuc.IZMAX = Z_count;
}

void initistest_();

void _hadgen_initis() {
    printf("OLOLOLO!!!\n");
   // filling the TABNUC database
   // Not a good solution, but works

   // Initialization
   isotope_count = 0;
   Z_count = 0;
   IAD1 = 0; IAD2 = 0;
   abunuc.L110 = 110;   // Not needed, but for clearance
   int i=0, j=0;
   for(i=0; i<110; i++) {
      abunuc.NUCADR[i][0] = 0.;
      abunuc.NUCADR[i][1] = 0.;
   }
   for(i=0; i<500; i++) abunuc.ABUNOR[i] = 0.;
   // Data filling (from original TABNUC.DAT)
#define HADGEN_TABNUC_DATA_INCLUDE
#include "tabnuc.h"

   // END OF DATABASE
   // Computing probabilities for sampling
   // REFER TO ORIGINAL CODE Hadgen.f
   
   // NOTICE: All indices (i,j,AD1,AD2) are FORTRAN indices
   for(j=1; j<=abunuc.IZMAX; j++) {
      int AD1=abunuc.NUCADR[j-1][0];  // FORTRAN indices !
      int AD2=abunuc.NUCADR[j-1][1];
      for(i=AD1; i<=AD2; i++) 
         abunuc.ABUNOR[i-1] = abunuc.PERC[i-1]; // FORTRAN indices !
      for(i=AD1; i<=AD2-1; i++)
         abunuc.ABUNOR[i] = abunuc.ABUNOR[i-1] + abunuc.ABUNOR[i];
      for(i=AD1; i<=AD2; i++)
         abunuc.ABUNOR[i-1]=abunuc.ABUNOR[i-1]/abunuc.ABUNOR[AD2-1];
   }
   // ALL DONE
//   initistest_(); // debug
}
