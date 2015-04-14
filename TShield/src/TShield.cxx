// @(#)root/transport/shield:$Id: TShield.cxx
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
#include "shield.h"     // link to $SIMPATH/generators/shield/shield.h
                        // includes all the lib Interface
#include "TShieldGeometry.h"

#include <cstdlib>
#include <time.h>

// class TShield implementation
TShield *TShield::fInstance = 0;


//______________________________________________________________________________
TShield::TShieldClean::~TShieldClean() {
   // delete the unique copy of TShield
   // if exists
   if (TShield::fInstance)
      delete TShield::fInstance;
}
//______________________________________________________________________________


ClassImp(TShield)

//______________________________________________________________________________
TShield::TShield()
{
    // Check if another instance of TShield already exists
    if (fInstance) {
        Fatal("TShield", "FATAL ERROR: Another instance of TShield class already exists");
    }
    fInstance = this;
    // create static cleaner
    static TShieldClean cleaner;
    
    SetAutoSeed();
    // Create arrays for particles storage
    fParticlesFlyOut = new TClonesArray("TParticle", ParticlesNumFlyOutMax);
    fParticlesAbsorbed = new TClonesArray("TParticle", ParticlesNumAbsorbedMax);
    fGeometry = 0;
    // init the random generator
    srand(time(0));
    // init callbacks for tree generation and geometry
    InitCallbacks();

//   fParticles = new TClonesArray("TParticle", 
//      shield_get_max_stp() + shield_get_max_snu()); 

    // Setting all default
    shield_set_defaults();
}

TShield* TShield::Instance() {
   return fInstance ? fInstance : new TShield();
}

//______________________________________________________________________________
TShield::~TShield() 
{
    // delete TClonesArray's
    delete fParticlesFlyOut;
    delete fParticlesAbsorbed;
    // if deleted, clear the instance
    TShield::fInstance = 0;
    // useful if needed to create a new instance of TShield in future
    // actually TShield object should be automatically cleaned
    // but user can also delete it manually
}

//______________________________________________________________________________
void TShield::GenerateEvent() {
    // Clean particles arrays
    fParticlesFlyOut->Clear();
    fParticlesAbsorbed->Clear();
    ParticlesNumFlyOut = 0;
    ParticlesNumAbsorbed = 0;
    
    shield_clean_geometry();
    Geant4ShieldData geometryData = convertTGeantToShield(fGeometry);
    for(int i = 0; i<geometryData.size();i++){
        int countBodies = geometryData.at(i).bodyVector.size();
        // at "23.2.4 Class template vector" says, that elements of a vector are stored contiguously.
        int* zones = &(geometryData.at(i).zoneVector[0]);
        SGeoBody* bodies = &(geometryData.at(i).bodyVector[0]); 
        if(shield_add_zone( countBodies, zones, bodies, geometryData.at(i).medium) == -1){
            printf("An error was accured");
        }
    }
    
    shield_run();
}

//______________________________________________________________________________
void TShield::SetGeometry(TGeoManager *geom) {
    fGeometry = geom;
}

//______________________________________________________________________________
void TShield::PrintGeometry() {
    tgeanttoshield::printGeant4ShieldData(convertTGeantToShield(fGeometry));
}

//______________________________________________________________________________
TGeoManager *TShield::GetGeometry() {
    return fGeometry;
}
