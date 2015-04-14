#ifndef CONVERT_GEANT4_TO_SHIELD_H
#define CONVERT_GEANT4_TO_SHIELD_H

#include <list>
#include <vector>
#include "shield.h"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "cmath"

#include "G4RunManager.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4CutTubs.hh"
#include "G4Cons.hh"
#include "G4Para.hh"
#include "G4Trd.hh"
#include "G4Trap.hh"
#include "G4Sphere.hh"
#include "G4Orb.hh"
#include "G4Torus.hh"
#include "G4Polycone.hh"
#include "G4Polyhedra.hh"
#include "G4EllipticalTube.hh"
#include "G4Ellipsoid.hh"
#include "G4Tet.hh"
#include "G4GenericTrap.hh"


// Output type (tuple) with data for adding to shield
struct Geant4ShieldElement {
    std::vector<int> zoneVector;
    std::vector<SGeoBody> bodyVector;
    MediumData medium;
};
typedef std::vector<Geant4ShieldElement> Geant4ShieldData;
Geant4ShieldData convertGeant4ToShield(G4RunManager *runManager);

namespace geant4toshield {
// Element for using in boolean calculations at geometry.
typedef std::pair<int, SGeoBody> zoneElement;
// vector zoneList is:
// zone, written as
// "(A1 A2 A3 A4) OR (B1 B2 B3 B4) OR (C1 C2 C3 C4)"
// where "(A1 A2)" is "A1 AND A2"
// at shield input files
// at this structure is:
// +---+-------------+
// | 1 | A1 A2 A3 A4 |
// +---+-------------+
// | 2 | B1 B2 B3 B4 |
// +---+-------------+
// | 3 | C1 C2 C3 C4 |
// +---+-------------+
typedef std::vector<std::vector<zoneElement> > zoneList;
// Type for one zone for Shield
typedef std::vector<std::pair<zoneList, MediumData> > zoneData;

// Boolean algebra function on zoneElement and zoneData
zoneElement notElement(zoneElement el);
zoneList notZone(zoneList list);
zoneList orZone(zoneList list1, zoneList list2);
zoneList andZone(zoneList list1, zoneList list2);

// Function for creating correct cone or cylinder with correct order of variables
SGeoBody createCone(G4ThreeVector startPoint, G4ThreeVector vectorToEnd, float rStart, float rEnd);
// Function for creating zoneList for bodies that sectioned along the z axis.
// Using with andZone functions for adding output value (that is substraction  zoneList value) to zone from full body.
zoneList cutByPhi(G4ThreeVector currTranslation, G4RotationMatrix currRotation, float sPhi, float dPhi, float halfX, float halfY, float halfZ);
// Functions for creating zones for primitive Geant4 bodies.
inline std::pair<zoneList, zoneList> boxToZones(G4Box *box, G4ThreeVector currTranslation, G4RotationMatrix currRotation, double scale);
inline std::pair<zoneList, zoneList> tubeToZones(G4Tubs *tube, G4ThreeVector currTranslation, G4RotationMatrix currRotation, double scale);
inline std::pair<zoneList, zoneList> coneToZones(G4Cons *cone, G4ThreeVector currTranslation, G4RotationMatrix currRotation, double scale);
inline std::pair<zoneList, zoneList> trdToZones(G4Trd *trd, G4ThreeVector currTranslation, G4RotationMatrix currRotation, double scale);
inline std::pair<zoneList, zoneList> sphereToZones(G4Sphere *sphere, G4ThreeVector currTranslation, G4RotationMatrix currRotation, double scale);
inline std::pair<zoneList, zoneList> orbToZones(G4Orb *orb, G4ThreeVector currTranslation, G4RotationMatrix currRotation, double scale);
// Converting polyhedra to shield geometry as substraction of external and internal polyhedras.
// Every plane presented as substraction of Shield's "wedge": One wedge is wedge from central line axis,
// other is wedge from one vertex with larger radius to wedge to smaller one.
inline std::pair<zoneList, zoneList> polyconeToZones(G4Polycone *polycone, G4ThreeVector currTranslation, G4RotationMatrix currRotation, double scale);
// Converting polyhedra to shield geometry as substraction of external and internal polyhedras.
// Every plane presented as substraction of Shield's "wedge": One wedge is wedge from central line axis,
// other is wedge from one vertex with larger radius to wedge to smaller one.
inline std::pair<zoneList, zoneList> polyhedraToZones(G4Polyhedra *polyhedra, G4ThreeVector currTranslation, G4RotationMatrix currRotation, double scale);
inline std::pair<zoneList, zoneList> ellipticalTubeToZones(G4EllipticalTube *tube, G4ThreeVector currTranslation, G4RotationMatrix currRotation, double scale);
inline std::pair<zoneList, zoneList> ellipsoidToZones(G4Ellipsoid *ell, G4ThreeVector currTranslation, G4RotationMatrix currRotation, double scale);
inline std::pair<zoneList, zoneList> genericTrapTubeToZones(G4GenericTrap *gtrap, G4ThreeVector currTranslation, G4RotationMatrix currRotation, double scale);

//Output value:
//  First pair: first is zone from body, second is outer shell of body (for substraction from mother volume)
//  Second pair: current rotation matrix and current translation vector for generating daughter volumes
inline double getDefaultScale(){return 1/cm;}
std::pair <  std::pair<zoneList, zoneList>, std::pair<G4RotationMatrix, G4ThreeVector> >
    getZoneFromBody(G4VPhysicalVolume *physBody, G4RotationMatrix oldRotation, G4ThreeVector oldTranslation, double scale = getDefaultScale());
std::pair<zoneList, MediumData> addOuterVacuum(G4VPhysicalVolume *world);
    MediumData getMedium(G4VPhysicalVolume *body);
std::pair < zoneList, std::pair<G4RotationMatrix, G4ThreeVector> >
convertToSGeoZone(G4VPhysicalVolume *physBody, G4RotationMatrix oldRotation, G4ThreeVector oldTranslation);
zoneData convertBodyRecursive(G4VPhysicalVolume *physBody, G4RotationMatrix oldRotation = G4RotationMatrix(), G4ThreeVector oldTranslation = G4ThreeVector());
Geant4ShieldData zoneDataToShield(zoneData data);

//Print used values
void printElement(Element element);
void printMediumData(MediumData medium);
void printSGeoBody(SGeoBody body);
void printSGeoZone(SGeoZone zone);
void printGeant4ShieldData (Geant4ShieldData data);

} //end of namespace geant4toshield

#endif
