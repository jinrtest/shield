#define SHIELD_LIB_INTERNAL
#include "shield.h"
#include "shield_common_blocks.h"
#include <stdio.h>

//TODO change "!=" to "notEqual(float a, float b) {return (fabs(a-b)<1E-7);}"
struct {
    int count;
    struct MediumData mediums[48];
} SMedium = {0, {}};

int _check_type(int type) {
    return ((type > 0) && (type < 4));
}

int _check_n_el(int n) {
    return ((n > 0) && (n <= 24));
}
int _isElementsEqual(struct Element element1, struct Element element2) {
    if (element1.Nuclid != element2.Nuclid ||
            element1.Conc != element2.Conc ||
            element1.Density != element2.Density ||
            element1.Z != element2.Z ||
            element1.A != element2.A ||
            element1.PureDensity != element2.PureDensity ||
            element1.ionEv != element2.ionEv) {
        return 0;
    }
    return 1;
}
int _isMediumEqual(struct MediumData medium1, struct MediumData medium2) {
    if (medium1.nType != medium2.nType ||
            medium1.nChemEl != medium2.nChemEl ||
            medium1.Rho != medium2.Rho) {
        return 0;
    }
    int i;
    for (i = 0; i < medium1.nChemEl; i++) {
        if (_isElementsEqual(medium1.Elements[i], medium2.Elements[i]) != 1) {
            return 0;
        }
    }
    return 1;
}

void shield_send_medium(struct MediumData medium, int MediaNum) {
    int j = MediaNum; // j is Fortran index
    macrin.NELEMD[j][0] = medium.nType;
    int NELEM = medium.nChemEl;
    if (medium.nType == 1) {
        NELEM = 1;
    }
    macrin.NELEMD[j][1] = NELEM;
    float RHO = medium.Rho;
    macrin.RHOMED[j] = RHO;

    int k;
    struct Element* el;
    for (k = 0; k < medium.nChemEl; k++) {
        el = &(medium.Elements[k]);
        int Nuclid = el->Nuclid;
        macrin.ELMIX[j][k][0] = (float)(el->Nuclid);
        macrin.ELMIX[j][k][1] = el->Conc;
        macrin.ELMIX[j][k][2] = el->Density;
        macrin.ELMIX[j][k][3] = (el->Z!=0) ? el->Z : infel.ZNUC[Nuclid - 1];  // was only infel.ZNUC[Nuclid - 1]. Why?
        macrin.ELMIX[j][k][4] = (el->A!=0) ? el->A : infel.ATWEI[Nuclid - 1]; // was only infel.ATWEI[Nuclid - 1]. Why?
        macrin.ELMIX[j][k][5] = (el->PureDensity != 0)? el->PureDensity : infel.DENS[Nuclid - 1]; // Fortran index!
        macrin.ELMIX[j][k][6] = (el->ionEv != 0) ? el->ionEv : infel.RIEV[Nuclid - 1];
        macrin.ELMIX[j][k][7] = el->reserved;
    }

    const float AVOG = 6.022169e-4;
    float RMOL;
    switch (medium.nType) {
        case 1:
            RHO = macrin.ELMIX[j][0][5];
            macrin.RHOMED[j] = RHO;
            macrin.ELMIX[j][0][1] = AVOG * (RHO / macrin.ELMIX[j][0][4]);
            macrin.ELMIX[j][0][2] = RHO;
            break;
        case 2: // accidentally medium 2 and 3 coincide
        case 3:
            RMOL = 0;
            for (k = 0; k < NELEM; k++) {
                RMOL += macrin.ELMIX[j][k][1] * macrin.ELMIX[j][k][4];
            }
            for (k = 0; k < NELEM; k++) {
                macrin.ELMIX[j][k][1] = (macrin.ELMIX[j][k][1] * RHO * AVOG) / RMOL;
                macrin.ELMIX[j][k][2] = (macrin.ELMIX[j][k][1] * macrin.ELMIX[j][k][4] * RHO) / RMOL;
            }
            break;
            /*case 3:
                float AVER = 0;
                for (k = 0; k < NELEM; k++) {
                    AVER += macrin.ELMIX[j][k][1] * macrin.ELMIX[j][k][4];
                }
                for (k = 0; k < NELEM; k++) {
                    macrin.ELMIX[j][k][1] = (macrin.ELMIX[j][k][1] * RHO * AVOG) / AVER;
                    macrin.ELMIX[j][k][2] = (macrin.ELMIX[j][k][1] * macrin.ELMIX[j][k][4] * RHO) / AVEr;
                }
                break;*/
        case 4: //Not completed in my version of shield
        default:
            break;
    }

    // success
    //MediaNum ++;
    //return MediaNum;
}

int shield_add_medium(struct MediumData medium) {
    int maxNum = 24; //Constant from fortran code
    int i;
    int mednum = -1;
    for (i = 0; i < SMedium.count; i++) {
        if (_isMediumEqual(medium, SMedium.mediums[i])) {
            mednum = i;
            break;
        }
    }
    if (mednum != -1) {
        return mednum;
    }
    if (SMedium.count >= maxNum) {
        printf("Maximum number of mediums. Only mediums #0 and #1000 now are allowed");
        return -1;
    } else if (SMedium.count == maxNum - 1) {
        printf("Maximum number of mediums are obtained. Only mediums #0 and #1000 fuhrer are allowed");
        return -1;
    }

    if (!_check_type(medium.nType)) {
        printf("Warning: Media type for medium %i is invalid (type = %i)\n\tMedium ignored.\n",medium.nType);
        return -1;
    }
    if (!_check_n_el(medium.nChemEl)) {
        printf("Warning: number of chem. el. for medium is invalid (= %i)\n\tMedium ignored.\n",medium.nChemEl);
        return -1;
    }

    SMedium.mediums[SMedium.count] = medium;
    return SMedium.count++;
}

void shield_init_medium(void) {
    int i;
    macrin.NUMMED = SMedium.count;
    for (i = 0; i < SMedium.count; i++) {
        shield_send_medium(SMedium.mediums[i], i);
    }
}
void shield_clean_medium(void) {
    SMedium.count=0;
}
