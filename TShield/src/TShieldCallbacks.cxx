// @(#)root/transport/shield:$Id: TShieldCallbacks.cxx
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
#include "TGeoManager.h"
#include "TGeoNode.h"
#include "TGeoVolume.h"
#include "TGeoMedium.h"

float TShield::GeoDistCallback(float x, float y, float z, 
                       float vx, float vy, float vz)
{   
    TShield *fShield = TShield::Instance();
    if (fShield->fGeometry != 0)
        return fShield->FindBoundary(x, y, z, vx, vy, vz);
    else {
        std :: cout << "TShield FATAL ERROR: TGeoManager not set" << std :: endl;
        return -1;
    }
}

int TShield::GeoNextCallback(float x, float y, float z, 
                     float vx, float vy, float vz)
{
    TShield *fShield = TShield::Instance();
    if (fShield->fGeometry != 0) 
        return fShield->GetNextVolume(x, y, z, vx, vy, vz);
    else {
        std :: cout << "TShield FATAL ERROR: TGeoManager not set" << std :: endl;
        return -1;
    }
}

void TShield::TreeCallback(shield_tree_node *node)
{
    switch(node->event) {
        case EVENT_FLYOUT:
            TShield::Instance()->AddFliedOutParticle(node);
            break;
        case EVENT_ABSORPTION:
            TShield::Instance()->AddAbsorbedParticle(node);
            break;
        default:
            break;
    }
    // anyway, all particles are always registered in the tree
    TShield::Instance()->AddTreeParticle(node);
}

void TShield::InitCallbacks()
{
    shield_set_tree_callback(&TShield::TreeCallback);
    shield_geo_set_dist_callback(&TShield::GeoDistCallback);
    shield_geo_set_next_callback(&TShield::GeoNextCallback);
}

void TShield::AddFliedOutParticle(shield_tree_node *node)
{
    if (ParticlesNumFlyOut < ParticlesNumFlyOutMax) {
        TClonesArray &a = *fParticlesFlyOut;
        new (a[ParticlesNumFlyOut]) TParticle(node->pdg, 0,0,0,0,0,
            node->px_z, node->py_z, node->pz_z, node->tz + node->weight,
            node->xz, node->yz, node->zz, 0);
        ParticlesNumFlyOut ++;
    } else {
        std::cout << "TShield ERROR: TClonesArray fParticlesFlyOut overfull" << std::endl;
    }
}

void TShield::AddAbsorbedParticle(shield_tree_node *node)
{
    if (ParticlesNumAbsorbed < ParticlesNumAbsorbedMax) {
        TClonesArray &a = *fParticlesAbsorbed;
        new (a[ParticlesNumAbsorbed]) TParticle(node->pdg, 0,0,0,0,0,
            node->px_z, node->py_z, node->pz_z, node->tz + node->weight,
            node->xz, node->yz, node->zz, 0);
        ParticlesNumAbsorbed ++;
    } else {
        std::cout << "TShield ERROR: TClonesArray fParticlesAbsorbed overfull" << std::endl;
    }
}

void TShield::AddTreeParticle(shield_tree_node *node)
{
    // TTree support should be added
}

float TShield::FindBoundary(float x, float y, float z, float vx, float vy, float vz)
{
    TGeoNode *node = fGeometry->FindNode(x, y, z);

}

int TShield::GetNextVolume(float x, float y, float z, float vx, float vy, float vz)
{
    TGeoNode *node = fGeometry->FindNode(x, y, z);    
}

