// @(#)root/transport/shield:$Id: TShieldOpt.cxx
// Authors: Alexander Timofeev 29/02/2012
//          Dmitry    Sosnov      02/2015

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// TShield                                                                    //
// is a Interface (wrapper) class for Shield transport code                   //
// Original code of SHIELD package                                            //
// is written by N.M.Sobolevsky                                               //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


#include "TShield.h"
#include "shield.h"
#include "hadgen.h"
#include <iostream>
#include <cmath>

void TShield::SetAutoSeed(bool flag)
{
    autoseed = flag;
}

void TShield::SetRandomSeed(Int_t seed)
{
    if (seed < 0) seed = - seed;
    shield_set_randomseed(seed);
}

Int_t TShield::GetRandomSeed()
{
    return shield_get_randomseed();
}

void TShield::SetLuxCount(Int_t luxcnt)
{
    shield_set_luxcount(luxcnt);
}

Int_t TShield::GetLuxCount()
{
    return shield_get_luxcount();
}

void TShield::SetLCASC(Int_t lcasc)
{
    shield_set_lcasc(lcasc);
}

Int_t TShield::GetLCASC()
{
    return shield_get_lcasc();
}

void TShield::SetIncidentParticle(Int_t jpart)
{
    shield_set_incidentparticle(jpart);
}

Int_t TShield::GetIncidentParticle()
{
    return shield_get_incidentparticle();
}

void TShield::SetEnergy(Float_t energy) 
{
    shield_set_energy(energy);
}

Float_t TShield::GetEnergy()
{
    return shield_get_energy();
}

void TShield::SetAProj(Float_t aproj)
{
    shield_set_aproj(aproj);
}

void TShield::SetZProj(Float_t zproj)
{
    shield_set_zproj(zproj);
}

Float_t TShield::GetAProj()
{
    return shield_get_aproj();
}

Float_t TShield::GetZProj()
{
    return shield_get_zproj();
}

void TShield::SetStartPoint(Float_t x, Float_t y, Float_t z)
{
    shield_set_runpoint(x, y, z);
}

void TShield::SetDirectionVector(Float_t cx, Float_t cy, Float_t cz)
{
    Float_t norm = sqrt(cx*cx + cy*cy + cz*cz);
    cx /= norm; cy /= norm; cz /= norm;
    Float_t sin_theta = sqrt(1 - cz*cz);
    if (sin_theta != 0.0) {
    shield_set_rundirection(cz, cy / sin_theta, cx / sin_theta);
    } else {
    shield_set_rundirection(1,0,1);
    }
}

void TShield::SetDirection(Float_t cos_theta, Float_t sin_phi, Float_t cos_phi)
{
    shield_set_rundirection(cos_theta, sin_phi, cos_phi);
}