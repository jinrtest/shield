#define SHIELD_LIB_INTERNAL
#include "shield.h"
#include "shield_common_blocks.h"

void shield_set_iso()
{
    isoset.NISTRU = 0;
    
#include "isotin.data"
}

void shield_add_iso(int iso1, int iso2, char tdec[10])
{
    int i = isoset.NISTRU;
    isoset.NISTRU ++;    
    isoset.ISOSET[i][0] = iso1;
    isoset.ISOSET[i][1] = iso2;
    int k;
    for (k=0; k<10; k++) isoset.TDEC[i][k] = tdec[k];
}
