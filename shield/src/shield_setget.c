#define SHIELD_LIB_INTERNAL
#include "shield.h"
#include "hadgen.h"
#include "shield_common_blocks.h"

void shield_set_randomseed(int seed)
{
    hadgen_set_randomseed(seed);
}

void shield_set_luxcount(int luxcnt)
{
    hadgen_set_luxcount(luxcnt);
}

void shield_set_incidentparticle(int jpart)
{
    t0jpr0.JPART0 = jpart;
}

void shield_set_energy(float energy)
{
    t0jpr0.TMAX0 = energy;
}

void shield_set_lcasc(int lcasc)
{
    debug.LCASC = lcasc;
}

void shield_set_statistics(int nstat)
{
    other.NSTAT = nstat;
}

void shield_set_nsave(int nsave)
{

}

void shield_set_aproj(float a)
{
    hipro0.APROJ = a;
}

void shield_set_zproj(float z)
{
    hipro0.ZPROJ = z;
}

int shield_get_randomseed()
{
    return hadgen_get_randomseed();
}

int shield_get_luxcount()
{
    return hadgen_get_luxcount();
}

int shield_get_incidentparticle()
{
    return t0jpr0.JPART0;
}

float shield_get_energy()
{
    return t0jpr0.TMAX0;
}

int shield_get_lcasc()
{
    return debug.LCASC;
}

int shield_get_statistics()
{
    return other.NSTAT;
}

int shield_get_nsave()
{
    return 1;
}

float shield_get_aproj()
{
    return hipro0.APROJ;
}

float shield_get_zproj()
{
    return hipro0.ZPROJ;
}

void shield_set_runpoint(float x, float y, float z)
{
    wrpr.WP[0] = x;
    wrpr.WP[1] = y;
    wrpr.WP[2] = z;    
}

void shield_set_rundirection(float c1, float s2, float c3)
{
    wrpr.WP[3] = c1;
    wrpr.WP[4] = s2;
    wrpr.WP[5] = c3;
}

