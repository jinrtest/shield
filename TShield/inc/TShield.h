// @(#)root/transport/shield:$Id: TShield.h
// Authors: Alexander Timofeev   29/02/2012
//          Dmitry    Sosnov        02/2015

#ifndef ROOT_Shield_header
#define ROOT_Shield_header

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TShield                                                              //
// is a wrapper class for SHIELD transport code                         //
// original code for SHIELD is                                          //
// written by N.M.Sobolevsky                                            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "TObjArray.h"
#include "TClonesArray.h"
#include "TGenerator.h"
#include "TParticle.h"
#include "TGeoManager.h"
#include "TTree.h"
#include "shield.h"
#include "stdio.h"

class TShield {
    ClassDef(TShield, 1)	// TShield
    
protected:
    // auto-cleanup mechanism
    static TShield *fInstance;
    class TShieldClean {
    public:
        TShieldClean() {}
        ~TShieldClean();
    };
    friend class TShieldClean;


    TClonesArray *fParticlesFlyOut, *fParticlesAbsorbed;
    Int_t ParticlesNumFlyOut, ParticlesNumAbsorbed;
    static const int ParticlesNumFlyOutMax = 10000;    // TODO
    static const int ParticlesNumAbsorbedMax = 10000;  // TODO
    bool autoseed;

    TGeoManager *fGeometry;
public:
    const TClonesArray *GetFlyOutArray() { return fParticlesFlyOut; }
    const TClonesArray *GetAbsorbedArray() { return fParticlesAbsorbed; }
    Int_t GetFlyOutNum() { return ParticlesNumFlyOut; }
    Int_t GetAbsorbedNum() { return ParticlesNumAbsorbed; }
private:	
   // TShield class should be never copied
    TShield(const TShield& s) {}
    TShield& operator=(const TShield&) { return *this; }
    
    void InitCallbacks();
    
    // static callbacks
    static int GeoNextCallback(float, float, float, float, float, float);
    static float GeoDistCallback(float, float, float, float, float, float);
    static void TreeCallback(shield_tree_node *);
    // non-static callbacks
    // to store the tree
    void AddTreeParticle(shield_tree_node *);
    // to store flied-out particles
    void AddFliedOutParticle(shield_tree_node *);
    // to store absorbed particles
    void AddAbsorbedParticle(shield_tree_node *);
    // to track geometry
    float FindBoundary(float x, float y, float z, float vx, float vy, float vz);
    int GetNextVolume(float x, float y, float z, float vx, float vy, float vz);    

    // Tree
    TTree *tree;

public:
    TShield();		// creates a unique instance of TShield
    virtual ~TShield();
	
    static TShield * Instance();

   // Random generator state
    void SetRandomSeed(Int_t seed);
    void SetRandomSeed();
    Int_t GetRandomSeed();

    void SetAutoSeed(bool flag = false);   // initially false    

    void SetGeometry(TGeoManager *geom);
    TGeoManager *GetGeometry();

protected:
   // the following functions should not be used
    void SetLuxCount(Int_t luxcnt);
    Int_t GetLuxCount();
    void SetLCASC(Int_t lcasc = 0);
    Int_t GetLCASC();

public:
   // Inreac properties
    void SetIncidentParticle(Int_t jpart);  // these functions may be more smart 
    Int_t GetIncidentParticle();


    //void SetNuclid(Int_t nuclid);
    //Int_t GetNuclid();


    void SetEnergy(Float_t energy);
    Float_t GetEnergy();
 
   // A and Z of projectile ion
   // notice that these params are set automaticaly for 
   // Deuteron, Tritium, He3 and alpha particles
    void SetAProj(Float_t aproj);
    Float_t GetAProj();
    void SetZProj(Float_t zproj);
    Float_t GetZProj();
    
    void SetStartPoint(Float_t x, Float_t y, Float_t z);
    void SetDirection(Float_t cos_theta, Float_t sin_phi, Float_t cos_phi);
    void SetDirectionVector(Float_t fx, Float_t cy, Float_t cz);
    void GenerateEvent();

    void PrintGeometry();

};


#endif
