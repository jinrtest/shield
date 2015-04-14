#define HADGEN_LIB_INTERNAL
#include "hadgen.h"
#include "hadgen_common_blocks.h"

static int Counter=0;

void _hdg_setdky(int IRES, int ITYPE, int IMODE1, int IMODE2, 
                 int IMODE3, int IMODE4, int IMODE5, double BR) {
   atabdata.IRES[Counter] = IRES;
   atabdata.IMODE[Counter][0] = IMODE1;
   atabdata.IMODE[Counter][1] = IMODE2;
   atabdata.IMODE[Counter][2] = IMODE3;
   atabdata.IMODE[Counter][3] = IMODE4;
   atabdata.IMODE[Counter][4] = IMODE5;
   atabdata.ITYPE[Counter] = ITYPE;
   atabdata.BR[Counter] = BR;
   Counter++;
}

void _hadgen_setdky() {
   // cleaning
   int i, j;
   for(i=0; i<600; i++) {
      for(j=0; j<5; j++) atabdata.IMODE[i][j] = 0;
      atabdata.IRES[i] = 0;
      atabdata.ITYPE[i] = 0;
      atabdata.BR[i] = 0.0;
   }
   Counter=0;
   // filling the database
#include "atab.h"
}

void dkytest_()
{
    
}