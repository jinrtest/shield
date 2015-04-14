// Converting geometry from TGeoManager's to SHIELD's.
// Authors: Dmitry    Sosnov      02/2015

#include "TShieldGeometry.h"

Geant4ShieldData convertTGeantToShield(TGeoManager* geoMan) {
    using namespace tgeanttoshield;
    TGeoNode* world = geoMan->GetTopNode();
    zoneData zoneWorld = convertBodyRecursive(world);
    zoneWorld.push_back(addOuterVacuum(world));
    return zoneDataToShield(zoneWorld);
}

bool operator<(const tgeanttoshield::zoneElement z1, const tgeanttoshield::zoneElement z2){
    if(z1.first!=z2.first){return (z1.first>z2.first);}
    if(z1.second.type!=z2.second.type){return (z1.second.type<z2.second.type);}
    for(int k = 0; k<36; ++k){
        if(z1.second.parameters[k]!=z2.second.parameters[k]){
            return (z1.second.parameters[k]<z2.second.parameters[k]);
        }
    }
    return false;
}

//Operators for comparsion floating point values
bool doubleNE(const double &l, const double &r){
    return (abs(l-r)>1E-15);
}
bool doubleEQ(const double &l, const double &r){
    return !doubleNE(l,r);
}
bool doubleLT(const double &l, const double &r){
    return (doubleNE(l,r) && (l<r));
}
bool doubleLE(const double &l, const double &r){
    return ((l<r) || doubleEQ(l,r));
}
bool doubleGT(const double &l, const double &r){
    return !doubleLE(l,r);
}
bool doubleGE(const double &l, const double &r){
    return !doubleLT(l,r);
}

namespace tgeanttoshield {
//-------------------------------------------------------------------------
zoneElement notElement(zoneElement el) {
    return std::pair<int, SGeoBody>(el.first * -1, el.second);
}
zoneList notZone(zoneList list) {
    int cur, count = 1;
    for (unsigned int i = 0; i < list.size(); ++i) {
        count *= list.at(i).size();
    }
    zoneList out;
    std::vector<zoneElement> tmp;
    for (int k = 0; k < count; ++k) {
        cur = k;
        tmp.clear();
        for (unsigned int i = 0; i < list.size(); ++i) {
            int cursize = list.at(i).size();
            tmp.push_back(notElement(list.at(i).at(cur % cursize)));
            cur /= cursize;
        }
        out.push_back(tmp);
    }
    return out;
}
zoneList orZone(zoneList list1, zoneList list2) {
    zoneList outList = list1;
    outList.insert(outList.end(), list2.begin(), list2.end());
    return outList;
}
zoneList andZone(zoneList list1, zoneList list2) {
    zoneList outList;
    std::vector<zoneElement> tmp;
    for (unsigned int i = 0; i < list1.size(); ++i) {
        for (unsigned int j = 0; j < list2.size(); ++j) {
            tmp = list1.at(i);
            tmp.insert(tmp.end(), list2.at(j).begin(), list2.at(j).end());
            outList.push_back(tmp);
        }
    }
    return outList;
}
//-------------------------------------------------------------------------
void printElement(Element element) {
    printf("Element: { %f,\t%f,\t%f,\t%f,\t%f,\t%f,\t%f}\n",
           element.Nuclid, element.Conc, element.Density,
           element.Z, element.A, element.PureDensity, element.ionEv);
}
void printMediumData(MediumData medium) {
    printf("Medium:\n");
    printf("\tType: %i\n", medium.nType);
    printf("\tRho: %f\n", medium.Rho);
    printf("\tCount of elements: %i\n", medium.nChemEl);
    for (int i = 0; i < medium.nChemEl; i++) {
        printf("\tElement %i: ", i);
        printElement(medium.Elements[i]);
    }
    printf("\n");
}
void printSGeoBody(SGeoBody body) {
    printf("Body: %i, { ", body.type);
    for (int i = 0; i < 36; i++) {
        printf("%f ", body.parameters[i]);
    }
    printf("}\n");
}
void printSGeoZone(SGeoZone zone) {
    printf("Zone:\n");
    printf("\tNumber of medium: %i\n", zone.mediaNumber);
    printf("\tCount of elements: %i\n", zone.countELements);
    printf("\t");
    for (int i = 0; i < zone.countELements; i++) {
        printf("%i ", zone.definition[2 * i]*zone.definition[2 * i + 1]);
    }
    printf("\n");
}
void printGeant4ShieldData (Geant4ShieldData data){
    for(unsigned int i = 0; i<data.size(); ++i){
        printf("Zones: ");
        for(unsigned int k = 0; k<data.at(i).zoneVector.size(); ++k){
            printf("%i ", data.at(i).zoneVector.at(k));
        }
        printf("\n  Bodies:\n");
        for(unsigned int k = 0; k<data.at(i).bodyVector.size(); ++k){
            printf("\t");
            printSGeoBody(data.at(i).bodyVector.at(k));
        }
        printf("  ");
        printMediumData(data.at(i).medium);
    }
}
//-------------------------------------------------------------------------

TGeoHMatrix operator*(TGeoHMatrix &matrixR, const TGeoMatrix &matrixL) {
    TGeoHMatrix mm = TGeoHMatrix(matrixR);
//     mm *= matrixL;
    mm.Multiply(&matrixL);
    return TGeoHMatrix(mm);
}
TGeoHMatrix operator*(TGeoHMatrix &matrixR, const TGeoHMatrix &matrixL) {
    TGeoHMatrix mm = TGeoHMatrix(matrixR);
    return mm *= matrixL;
}
TGeoHMatrix operator*(const TGeoMatrix &matrixR, const TGeoMatrix &matrixL) {
    TGeoHMatrix hmr = TGeoHMatrix(matrixR);
    return hmr * matrixL;
}
TGeoHMatrix operator*(const TGeoMatrix &matrixR, const TGeoHMatrix &matrixL) {
    TGeoHMatrix hmr = TGeoHMatrix(matrixR);
    return hmr * matrixL;
}
TGeoTranslation operator+(const TGeoTranslation &matrixR, const TGeoTranslation &matrixL) {
    TGeoTranslation m = TGeoTranslation(matrixR);
    m.Add(&matrixL);
    return TGeoTranslation(m);
}
TGeoTranslation operator-(const TGeoTranslation matrixR, const TGeoTranslation &matrixL) {
    return matrixR + matrixL.Inverse();
}
TGeoTranslation operator+(const TGeoTranslation matrixR, const TGeoMatrix &matrixL) {
    return matrixR + TGeoTranslation(matrixL);
}
TGeoTranslation operator*(const TGeoTranslation matrixR, const Double_t &scale) {
    Double_t* tmp1 = (Double_t*)matrixR.GetTranslation();
    return TGeoTranslation(tmp1[0]*scale,tmp1[1]*scale,tmp1[2]*scale);
}
void addVectorToElement(SGeoBody &body, unsigned int position, TGeoTranslation vector){
    body.parameters[position] = (float)(vector.GetTranslation()[0]);
    body.parameters[position+1] = (float)(vector.GetTranslation()[1]);
    body.parameters[position+2] = (float)(vector.GetTranslation()[2]);
}

//-------------------------------------------------------------------------
Geant4ShieldData zoneDataToShield(zoneData data) {
    Geant4ShieldData out;
    Geant4ShieldElement tmpElement;
    std::pair<zoneList, MediumData> currentPair;
    zoneList currentList;
    std::vector<SGeoBody> bodyVector;
    std::vector<int> zoneVector;
    std::set <zoneElement> curSet; //Set have only unique keys.
    for (unsigned int kk = 0; kk < data.size(); ++kk) {
        bodyVector.clear();
        zoneVector.clear();
        currentPair = data.at(kk);
        currentList = currentPair.first;
        for (unsigned int k = 0; k < currentList.size(); ++k) {
            curSet = std::set<zoneElement>(currentList.at(k).begin(), currentList.at(k).end()); //As fact, removing duplicates of zoneElements
            for (std::set<zoneElement>::const_iterator i = curSet.begin(); i != curSet.end(); ++i) {
                zoneVector.push_back(i->first);
                bodyVector.push_back(i->second);
            }
//             for (unsigned int i = 0; i < currentList.at(k).size(); ++i) { //Здесь только те тела, которые +1 или -1.
//                 zoneVector.push_back(currentList.at(k).at(i).first);
//                 bodyVector.push_back(currentList.at(k).at(i).second);
//             }
            zoneVector.push_back(0);
        }
        zoneVector.pop_back();
        tmpElement.medium = currentPair.second;
        tmpElement.zoneVector = zoneVector;
        tmpElement.bodyVector = bodyVector;
        out.push_back(tmpElement);
    }
    return out;
}
//-------------------------------------------------------------------------
zoneData convertBodyRecursive(TGeoNode *node) {
    return convertBodyRecursive(node, TGeoIdentity()).first;
}
std::pair<zoneData, zoneList> convertBodyRecursive(TGeoNode *node, TGeoHMatrix oldTransformation) {
    zoneData out;
    std::pair<std::pair<zoneList, zoneList>, TGeoHMatrix > tmp = getZoneFromBody(node, oldTransformation);
    TGeoHMatrix currMatrix = tmp.second;
    zoneList internalVolume = tmp.first.first;
    zoneList outerVolume = tmp.first.second;

    TGeoVolume* volume = node->GetVolume();
    int countOfDaughters = volume->GetNdaughters();
    TGeoNode *daughter;
    std::pair<zoneData, zoneList > daughterZones;
    for (int i = 0; i < countOfDaughters; ++i) {
        daughter = volume->GetNode(i);
        daughterZones = convertBodyRecursive(daughter, currMatrix);
        out.insert(out.end(), daughterZones.first.begin(), daughterZones.first.end()); //Insert daughter zones to output zones
        internalVolume = andZone(internalVolume, notZone(daughterZones.second)); //Only outel shell of body.
    }
//     out.push_back(std::make_pair(internalVolume, getMedium(node)));
    out.insert(out.begin(), std::make_pair(internalVolume, getMedium(node))); //Add internalVolume to first place at output.
    return std::make_pair(out,outerVolume);
}

std::pair<std::pair<zoneList, zoneList>, TGeoHMatrix > getZoneFromBody(TGeoNode *node, TGeoHMatrix oldTransformation, double scale) {
    TGeoHMatrix newTransformation = *(node->GetMatrix());
    TGeoHMatrix currTransformation = oldTransformation * newTransformation;
    TGeoShape *shape = node->GetVolume()->GetShape();
    return std::make_pair(getZoneFromShape(shape,currTransformation,scale),currTransformation);
}
    
std::pair<zoneList, zoneList> getZoneFromShape(TGeoShape *shape,TGeoHMatrix currTransformation, double scale) {
    std::pair<zoneList, zoneList> zonePair;
    if (dynamic_cast<TGeoTubeSeg*>(shape) != NULL) {
        zonePair = tubeSegToZones((TGeoTubeSeg*)shape, currTransformation, scale);
    } else if (dynamic_cast<TGeoTube*>(shape) != NULL) {
        zonePair = tubeToZones((TGeoTube*)shape, currTransformation, scale);
    } else if (dynamic_cast<TGeoCtub*>(shape) != NULL) {//TODO
        printf("TGeoCtub not implemented yet\n");
    } else if (dynamic_cast<TGeoConeSeg*>(shape) != NULL) {
        zonePair = coneSegToZones((TGeoConeSeg*)shape, currTransformation, scale);
    } else if (dynamic_cast<TGeoCone*>(shape) != NULL) {
        zonePair = coneToZones((TGeoCone*)shape, currTransformation, scale);
    } else if (dynamic_cast<TGeoPara*>(shape) != NULL) {//TODO
        printf("TGeoPara not implemented yet\n");
    } else if (dynamic_cast<TGeoTrd1*>(shape) != NULL) {
        zonePair = trd1ToZones((TGeoTrd1*)shape, currTransformation, scale);
    } else if (dynamic_cast<TGeoTrd2*>(shape) != NULL) {
        zonePair = trd2ToZones((TGeoTrd2*)shape, currTransformation, scale);
    } else if (dynamic_cast<TGeoTrap*>(shape) != NULL) { //TODO //I cannot understand angles and parameters from getters
        printf("TGeoTrap not implemented yet\n");
//     } else if (dynamic_cast<G4Orb *>(solidBody) != NULL) {
//         zonePair = orbToZones((G4Orb *)solidBody, currTranslation, currRotation, scale);
    } else if (dynamic_cast<TGeoSphere*>(shape) != NULL) {
        zonePair = sphereToZones((TGeoSphere*)shape, currTransformation, scale);
    } else if (dynamic_cast<TGeoPgon*>(shape) != NULL) {
        zonePair = polyhedraToZones((TGeoPgon*)shape, currTransformation, scale);
    } else if (dynamic_cast<TGeoPcon*>(shape) != NULL) {
        zonePair = polyconeToZones((TGeoPcon*)shape, currTransformation, scale);
    } else if (dynamic_cast<TGeoEltu*>(shape) != NULL) {
        zonePair = ellipticalTubeToZones((TGeoEltu*)shape, currTransformation, scale);
//     } else if (dynamic_cast<G4Ellipsoid*>(solidBody) != NULL) {
//         zonePair = ellipsoidToZones((G4Ellipsoid*)solidBody, currTranslation, currRotation, scale);
    } else if (dynamic_cast<TGeoArb8*>(shape) != NULL) {
        zonePair = arb8ToZones((TGeoArb8*)shape, currTransformation, scale);
//     } else if (dynamic_cast<G4Tet*>(solidBody) != NULL) { //TODO
//         printf("G4Tet not implemented yet\n");
    } else if (dynamic_cast<TGeoGtra*>(shape) != NULL) { //TODO
        printf("TGeoGtra not implemented yet\n");
    } else if (dynamic_cast<TGeoHype*>(shape) != NULL) { //TODO
        printf("TGeoHype not implemented yet\n");
    } else if (dynamic_cast<TGeoTorus*>(shape) != NULL) { //TODO
        printf("TGeoTorus not implemented yet\n");
    } else if (dynamic_cast<TGeoParaboloid*>(shape) != NULL) { //TODO
        printf("TGeoParaboloid not implemented yet\n");
    } else if (dynamic_cast<TGeoXtru*>(shape) != NULL) { //TODO
        printf("TGeoXtru not implemented yet\n");
    } else if (dynamic_cast<TGeoHalfSpace*>(shape) != NULL) { //TODO
        printf("TGeoHalfSpace not implemented yet\n");
    } else if (dynamic_cast<TGeoCompositeShape*>(shape) != NULL) { //TODO
        printf("TGeoCompositeShape not implemented yet\n");
    } else if (dynamic_cast<TGeoBBox*>(shape) != NULL) { //All of this shapes in subclasses of TGeoBBox
        zonePair = boxToZones((TGeoBBox*)shape, currTransformation, scale);
    } else {
        printf("ELSE");
    }
    return zonePair;
}

SGeoBody createCone(TGeoTranslation startPoint, TGeoTranslation vectorToEnd, float rStart, float rEnd) {
    SGeoBody outBody;
    if (doubleEQ(rStart,rEnd)) {
        outBody = {5, {0,0,0,0,0,0, rStart}};
        addVectorToElement(outBody,0,startPoint);
        addVectorToElement(outBody,3,vectorToEnd);
    } else {
        startPoint = (rStart > rEnd) ? startPoint : startPoint + vectorToEnd; //WARNING
        vectorToEnd = (rStart > rEnd) ? vectorToEnd : vectorToEnd.Inverse(); //WARNING
        float rMax = (rStart > rEnd) ? rStart : rEnd;
        float rMin = (rStart > rEnd) ? rEnd : rStart;
        outBody = {7, {0,0,0,0,0,0, rMax,rMin}};
        addVectorToElement(outBody,0,startPoint);
        addVectorToElement(outBody,3,vectorToEnd);
    }
    return outBody;
}
zoneList cutByPhi(TGeoHMatrix currTransformation,
                  float sPhi, float dPhi, float halfX, float halfY, float halfZ) {
    zoneList outList;
    TGeoHMatrix currRotation = TGeoHMatrix(currTransformation); currRotation.SetTranslation(kNullVector);
    TGeoTranslation currTranslation = TGeoTranslation(currTransformation);
    float ePhi = sPhi + dPhi;
    if(doubleGE(dPhi,360)) return {{}};
    TGeoHMatrix tmp = TGeoIdentity();tmp.RotateZ(sPhi);
    TGeoHMatrix innerRot = tmp * currRotation; //WARNING
    TGeoTranslation startFirstBox = currTranslation + TGeoTranslation(innerRot * TGeoTranslation(0, -halfY, -halfZ)); //WARNING
    TGeoTranslation vec11 = TGeoTranslation(innerRot * TGeoTranslation(-halfX, 0, 0));
    TGeoTranslation vec12 = TGeoTranslation(innerRot * TGeoTranslation(0, 2 * halfY, 0));
    TGeoTranslation vec13 = TGeoTranslation(innerRot * TGeoTranslation(0, 0, 2 * halfZ));
    SGeoBody box1 = {3, {}};
    addVectorToElement(box1,0,startFirstBox);
    addVectorToElement(box1,3,vec11);
    addVectorToElement(box1,6,vec12);
    addVectorToElement(box1,9,vec13);
    tmp = TGeoIdentity();tmp.RotateZ(ePhi);
    TGeoHMatrix outerRot = tmp * currRotation; //WARNING
    TGeoTranslation startSecondBox = currTranslation + TGeoTranslation(outerRot * TGeoTranslation(0, -halfY, -halfZ)); //WARNING
    TGeoTranslation vec21 = TGeoTranslation(outerRot * TGeoTranslation(+halfX, 0, 0));
    TGeoTranslation vec22 = TGeoTranslation(outerRot * TGeoTranslation(0, 2 * halfY, 0));
    TGeoTranslation vec23 = TGeoTranslation(outerRot * TGeoTranslation(0, 0, 2 * halfZ));
    SGeoBody box2 = {3, {}};
    addVectorToElement(box2,0,startSecondBox);
    addVectorToElement(box2,3,vec21);
    addVectorToElement(box2,6,vec22);
    addVectorToElement(box2,9,vec23);
    if (doubleLE(dPhi,180)) {
        outList = {{std::make_pair(-1, box1), std::make_pair(-1, box2)}};
    } else {
        outList = {{std::make_pair(-1, box1)}, {std::make_pair(-1, box2)}};
    }
    return outList;
}
inline std::pair<zoneList, zoneList> boxToZones(TGeoBBox *box,TGeoHMatrix currTransformation, double scale) {
    zoneList out, outerShell; //outer shell for conjunction at generating mother volume
    TGeoHMatrix currRotation = TGeoHMatrix(currTransformation); currRotation.SetTranslation(kNullVector);
    float x = (float) box->GetDX() * scale;
    float y = (float) box->GetDY() * scale;
    float z = (float) box->GetDZ() * scale;
    TGeoTranslation startVector = TGeoTranslation(currTransformation * TGeoTranslation(-x, -y, -z));
    TGeoTranslation vec1 = TGeoTranslation(currRotation * TGeoTranslation(2 * x, 0, 0));
    TGeoTranslation vec2 = TGeoTranslation(currRotation * TGeoTranslation(0, 2 * y, 0));
    TGeoTranslation vec3 = TGeoTranslation(currRotation * TGeoTranslation(0, 0, 2 * z));
    SGeoBody tmp = {3, {}};
    addVectorToElement(tmp,0,startVector);
    addVectorToElement(tmp,3,vec1);
    addVectorToElement(tmp,6,vec2);
    addVectorToElement(tmp,9,vec3);
    out = {{std::make_pair(1, tmp)}};
    outerShell = out;
    return std::make_pair(out, outerShell);
}
inline std::pair<zoneList, zoneList> tubeSegToZones(TGeoTubeSeg *tube,TGeoHMatrix currTransformation, double scale) {
    zoneList out, outerShell; //outer shell for conjunction at generating mother volume
    std::pair<zoneList, zoneList> tubeOut = tubeToZones(tube, currTransformation, scale);
    float sPhi = (float) tube->GetPhi1(); // As I understand, in gedrees
    float dPhi = (float) tube->GetPhi2() - sPhi;
    float rOut = (float) tube->GetRmax() * scale;
    float z = (float) tube->GetDz() * scale;
    zoneList phiZone = cutByPhi(currTransformation, sPhi, dPhi, rOut, rOut, z);
    out = andZone(tubeOut.first, phiZone);
    outerShell = andZone(tubeOut.second, phiZone);
    return std::make_pair(out, outerShell);
}
inline std::pair<zoneList, zoneList> tubeToZones(TGeoTube *tube,TGeoHMatrix currTransformation, double scale) {
    zoneList out, outerShell; //outer shell for conjunction at generating mother volume
    TGeoHMatrix currRotation = TGeoHMatrix(currTransformation); currRotation.SetTranslation(kNullVector);
    float rIn = (float) tube->GetRmin() * scale;
    float rOut = (float) tube->GetRmax() * scale;
    float z = (float) tube->GetDz() * scale;
    TGeoTranslation startTube = TGeoTranslation(currTransformation * TGeoTranslation(0, 0,-z));
    TGeoTranslation endTubeVec = TGeoTranslation(currRotation * TGeoTranslation(0, 0, 2 * z));
    SGeoBody innerTube = {5, {}};
    addVectorToElement(innerTube,0,startTube);
    addVectorToElement(innerTube,3,endTubeVec);
    innerTube.parameters[6] = rIn;
    SGeoBody outerTube = {5, {}};
    addVectorToElement(outerTube,0,startTube);
    addVectorToElement(outerTube,3,endTubeVec);
    outerTube.parameters[6] = rOut;
    out = {{std::make_pair(1, outerTube)}};
    outerShell = out;
    if (doubleNE(rIn,0)) {
        out = andZone(out, {{std::make_pair(-1, innerTube)}});
    }
    return std::make_pair(out, outerShell);
}
inline std::pair<zoneList, zoneList> coneSegToZones(TGeoConeSeg *cone,TGeoHMatrix currTransformation, double scale) {
    zoneList out, outerShell; //outer shell for conjunction at generating mother volume
    float rMax1 = (float)cone->GetRmax1() * scale;
    float rMax2 = (float)cone->GetRmax2() * scale;
    float z = (float)cone->GetDz() * scale;
    float sPhi = (float) cone->GetPhi1(); // As I understand, in gedrees
    float dPhi = (float) cone->GetPhi2() - sPhi;
    
    std::pair<zoneList, zoneList> coneOut = coneToZones(cone, currTransformation,scale);

    float rOut = (rMax1 > rMax2) ? rMax1 : rMax2;
    zoneList phiZone = cutByPhi(currTransformation, sPhi, dPhi, rOut, rOut, z);
    out = andZone(coneOut.first, phiZone);
    outerShell = andZone(coneOut.second, phiZone);
    return std::make_pair(out, outerShell);
}
inline std::pair<zoneList, zoneList> coneToZones(TGeoCone *cone,TGeoHMatrix currTransformation, double scale) {
    zoneList out, outerShell; //outer shell for conjunction at generating mother volume
    TGeoHMatrix currRotation = TGeoHMatrix(currTransformation); currRotation.SetTranslation(kNullVector);
    float rMin1 = (float)cone->GetRmin1() * scale;
    float rMin2 = (float)cone->GetRmin2() * scale;
    float rMax1 = (float)cone->GetRmax1() * scale;
    float rMax2 = (float)cone->GetRmax2() * scale;
    float z = (float)cone->GetDz() * scale;
    
    TGeoTranslation startTubeVector = TGeoTranslation(currTransformation * TGeoTranslation(0, 0, -z));
    TGeoTranslation endTubeVec = TGeoTranslation(currRotation * TGeoTranslation(0, 0, 2 * z));
    SGeoBody outerCone = createCone(startTubeVector, endTubeVec, rMax1, rMax2);
    SGeoBody innerCone = createCone(startTubeVector, endTubeVec, rMin1, rMin2);

    out = {{std::make_pair(1, outerCone)}};
    outerShell = out;
    if (doubleNE(rMin1,0) || doubleNE(rMin2,0)) {
        out = andZone(out, {{std::make_pair(-1, innerCone)}});
    }

    return std::make_pair(out, outerShell);
}
inline std::pair<zoneList, zoneList> trd1ToZones(TGeoTrd1 *trd,TGeoHMatrix currTransformation, double scale) {
    zoneList out, outerShell; //outer shell for conjunction at generating mother volume
    float dz = (float) trd->GetDz() * scale;
    float dx1 = (float) trd->GetDx1() * scale;
    float dy1 = (float) trd->GetDy() * scale;
    float dx2 = (float) trd->GetDx2() * scale;
    float dy2 = (float) trd->GetDy() * scale;
    TGeoTranslation v1 = TGeoTranslation(currTransformation * TGeoTranslation(-dx1, -dy1, -dz));
    TGeoTranslation v2 = TGeoTranslation(currTransformation * TGeoTranslation(-dx1, dy1, -dz));
    TGeoTranslation v3 = TGeoTranslation(currTransformation * TGeoTranslation(dx1, dy1, -dz));
    TGeoTranslation v4 = TGeoTranslation(currTransformation * TGeoTranslation(dx1, -dy1, -dz));
    TGeoTranslation v5 = TGeoTranslation(currTransformation * TGeoTranslation(-dx2, -dy2, dz));
    TGeoTranslation v6 = TGeoTranslation(currTransformation * TGeoTranslation(-dx2, dy2, dz));
    TGeoTranslation v7 = TGeoTranslation(currTransformation * TGeoTranslation(dx2, dy2, dz));
    TGeoTranslation v8 = TGeoTranslation(currTransformation * TGeoTranslation(dx2, -dy2, dz));
    SGeoBody tmp = {2, {}};
    addVectorToElement(tmp,0,v1);
    addVectorToElement(tmp,3,v2);
    addVectorToElement(tmp,6,v3);
    addVectorToElement(tmp,9,v4);
    addVectorToElement(tmp,12,v5);
    addVectorToElement(tmp,15,v6);
    addVectorToElement(tmp,18,v7);
    addVectorToElement(tmp,21,v8);
    out = {{std::make_pair(1, tmp)}};
    outerShell = out;
    return std::make_pair(out, outerShell);
}
inline std::pair<zoneList, zoneList> trd2ToZones(TGeoTrd2 *trd,TGeoHMatrix currTransformation, double scale) {
    zoneList out, outerShell; //outer shell for conjunction at generating mother volume
    float dz = (float) trd->GetDz() * scale;
    float dx1 = (float) trd->GetDx1() * scale;
    float dy1 = (float) trd->GetDy1() * scale;
    float dx2 = (float) trd->GetDx2() * scale;
    float dy2 = (float) trd->GetDy2() * scale;
    TGeoTranslation v1 = TGeoTranslation(currTransformation * TGeoTranslation(-dx1, -dy1, -dz));
    TGeoTranslation v2 = TGeoTranslation(currTransformation * TGeoTranslation(-dx1, dy1, -dz));
    TGeoTranslation v3 = TGeoTranslation(currTransformation * TGeoTranslation(dx1, dy1, -dz));
    TGeoTranslation v4 = TGeoTranslation(currTransformation * TGeoTranslation(dx1, -dy1, -dz));
    TGeoTranslation v5 = TGeoTranslation(currTransformation * TGeoTranslation(-dx2, -dy2, dz));
    TGeoTranslation v6 = TGeoTranslation(currTransformation * TGeoTranslation(-dx2, dy2, dz));
    TGeoTranslation v7 = TGeoTranslation(currTransformation * TGeoTranslation(dx2, dy2, dz));
    TGeoTranslation v8 = TGeoTranslation(currTransformation * TGeoTranslation(dx2, -dy2, dz));
    SGeoBody tmp = {2, {}};
    addVectorToElement(tmp,0,v1);
    addVectorToElement(tmp,3,v2);
    addVectorToElement(tmp,6,v3);
    addVectorToElement(tmp,9,v4);
    addVectorToElement(tmp,12,v5);
    addVectorToElement(tmp,15,v6);
    addVectorToElement(tmp,18,v7);
    addVectorToElement(tmp,21,v8);
    out = {{std::make_pair(1, tmp)}};
    outerShell = out;
    return std::make_pair(out, outerShell);
}
inline std::pair<zoneList, zoneList> sphereToZones(TGeoSphere *sphere,TGeoHMatrix currTransformation, double scale) {
    zoneList out, outerShell; //outer shell for conjunction at generating mother volume
    TGeoHMatrix currRotation = TGeoHMatrix(currTransformation); currRotation.SetTranslation(kNullVector);
    
    float rIn = (float) sphere->GetRmin() * scale;
    float rOut = (float) sphere->GetRmax() * scale;
    float sPhi = (float) sphere->GetPhi1();
    float dPhi = (float) sphere->GetPhi2() - sPhi;
    float sTheta = (float) sphere->GetTheta1() / (180/M_PI); //convert to radians becouse cmath functions work with them
    float dTheta = (float) sphere->GetTheta2() / (180/M_PI) - sTheta;
    float eTheta = sTheta + dTheta;
    SGeoBody innerSphere = {0, {}};
    innerSphere.parameters[3] = rIn;
    SGeoBody outerSphere = {0, {}};
    outerSphere.parameters[3] = rOut;
    {
        TGeoTranslation currTranslation = TGeoTranslation(currTransformation);
        addVectorToElement(innerSphere,0,currTranslation);
        addVectorToElement(outerSphere,0,currTranslation);
    }
    out = {{std::make_pair(1, outerSphere)}};
    outerShell = out;
    if (doubleNE(rIn,0)) {
        out = andZone(out, {{std::make_pair(-1, innerSphere)}});
    }
    zoneList phiZone = cutByPhi(currTransformation, sPhi, dPhi, rOut, rOut, rOut);
    out = andZone(out, phiZone);
    outerShell = andZone(outerShell, phiZone);

    if (doubleNE(sTheta,0) && doubleNE(sTheta,M_PI) && doubleNE(sTheta,M_PI/2.0)) {
        TGeoTranslation trc1Start = TGeoTranslation(currTransformation * TGeoTranslation(0, 0, -rOut * cos(sTheta)));
        TGeoTranslation trc1Vec = TGeoTranslation(currRotation * TGeoTranslation(0, 0, (rOut - rIn) * cos(sTheta)));
        SGeoBody trc1 = {7, {}};
        addVectorToElement(trc1,0,trc1Start);
        addVectorToElement(trc1,3,trc1Start);
        trc1.parameters[6] = (float)(rOut * tan(sTheta));
        trc1.parameters[7] = (float)(rIn * sin(sTheta));
        int mul = (sTheta < (M_PI/2.0)) ? -1 : 1;
        out = andZone(out, {{std::make_pair(mul, trc1)}});
        outerShell = andZone(outerShell, {{std::make_pair(mul, trc1)}});
    } else if (sTheta == (M_PI/2.0)) {
        TGeoTranslation db1s = TGeoTranslation(currTransformation * TGeoTranslation(-rOut, -rOut, 0));
        TGeoTranslation db1v1 = TGeoTranslation(currRotation * TGeoTranslation(2 * rOut, 0, 0));
        TGeoTranslation db1v2 = TGeoTranslation(currRotation * TGeoTranslation(0, 2 * rOut, 0));
        TGeoTranslation db1v3 = TGeoTranslation(currRotation * TGeoTranslation(0, 0, -rOut));
        SGeoBody db1 = {3, {}};
        addVectorToElement(db1,0,db1s);
        addVectorToElement(db1,3,db1v1);
        addVectorToElement(db1,6,db1v2);
        addVectorToElement(db1,9,db1v3);
        out = andZone(out, {{std::make_pair(-1, db1)}});
        outerShell = andZone(outerShell, {{std::make_pair(-1, db1)}});
    }

    if (doubleNE(eTheta,0) && doubleNE(eTheta,M_PI) && doubleNE(eTheta,M_PI/2.0)) {
        TGeoTranslation trc2Start = TGeoTranslation(currTransformation * TGeoTranslation(0, 0, -rOut * cos(eTheta)));
        TGeoTranslation trc2Vec = TGeoTranslation(currRotation * TGeoTranslation(0, 0, (rOut - rIn) * cos(eTheta)));
        SGeoBody trc2 = {7, {}};
        addVectorToElement(trc2,0,trc2Start);
        addVectorToElement(trc2,3,trc2Start);
        trc2.parameters[6] = (float)(rOut * tan(eTheta));
        trc2.parameters[7] = (float)(rIn * sin(eTheta));
        int mul = (eTheta > (M_PI/2.0)) ? -1 : 1;
        out = andZone(out, {{std::make_pair(mul, trc2)}});
        outerShell = andZone(outerShell, {{std::make_pair(mul, trc2)}});
    } else if (eTheta == (M_PI/2.0)) {
        TGeoTranslation db2s = TGeoTranslation(currTransformation * TGeoTranslation(-rOut, -rOut, 0));
        TGeoTranslation db2v1 = TGeoTranslation(currRotation * TGeoTranslation(2 * rOut, 0, 0));
        TGeoTranslation db2v2 = TGeoTranslation(currRotation * TGeoTranslation(0, 2 * rOut, 0));
        TGeoTranslation db2v3 = TGeoTranslation(currRotation * TGeoTranslation(0, 0, rOut));
        SGeoBody db2 = {3, {}};
        addVectorToElement(db2,0,db2s);
        addVectorToElement(db2,3,db2v1);
        addVectorToElement(db2,6,db2v2);
        addVectorToElement(db2,9,db2v3);
        out = andZone(out, {{std::make_pair(-1, db2)}});
        outerShell = andZone(outerShell, {{std::make_pair(-1, db2)}});
    }
    return std::make_pair(out, outerShell);
}
inline std::pair<zoneList, zoneList> polyconeToZones(TGeoPcon *polycone,TGeoHMatrix currTransformation, double scale) {
    zoneList out, outerShell; //outer shell for conjunction at generating mother volume
    TGeoHMatrix currRotation = TGeoHMatrix(currTransformation); currRotation.SetTranslation(kNullVector);
    float sPhi = (float)polycone->GetPhi1(); //In Pcon angle converted to degrees, and TMatrixes work with degrees. (?)
    float dPhi = (float)polycone->GetDphi();
    int nz = polycone->GetNz();
    float zCurrent, zPrev, rMaxCurrent, rMaxPrev, rMinCurrent, rMinPrev;
    double *zValues = polycone->GetZ();
    double *rMinValues = polycone->GetRmin();
    double *rMaxValues = polycone->GetRmax();
    float zMin = zValues[0] * scale, zMax = zValues[nz - 1] * scale, rMaxMax = rMaxValues[0] * scale;
    TGeoTranslation startPoint, endVec;
    zoneElement currentOuterCone, currentInnerCone;
    float zCenter = (zMax + zMin) / 2.0;
    for (int i = 1; i < nz; ++i) {
        zPrev = zCenter + (zValues[i - 1] - zCenter) * scale;
        zCurrent = zCenter + (zValues[i] - zCenter) * scale;
        if(zPrev == zCurrent)continue;
        rMaxPrev = rMaxValues[i - 1] * scale;
        rMaxCurrent = rMaxValues[i] * scale;
        rMinPrev = rMinValues[i - 1] * scale;
        rMinCurrent = rMinValues[i] * scale;
        startPoint = TGeoTranslation(currTransformation * TGeoTranslation(0, 0, zPrev));
        endVec = TGeoTranslation(currRotation * TGeoTranslation(0, 0, zCurrent - zPrev));
        rMaxMax = (rMaxCurrent > rMaxMax) ? rMaxCurrent : rMaxMax;
        if (doubleNE(rMaxPrev,0) || doubleNE(rMaxCurrent,0)) {
            currentOuterCone = std::make_pair(1, createCone(startPoint, endVec, rMaxPrev, rMaxCurrent));
            outerShell = orZone(outerShell, {{currentOuterCone}});
        }
        if (doubleNE(rMinPrev,0) || doubleNE(rMinCurrent,0)) {
            currentInnerCone = std::make_pair(1, createCone(startPoint, endVec, rMinPrev, rMinCurrent));
            out = orZone(out, {{currentInnerCone}});
        }
    }
    out = andZone(outerShell, notZone(out));
    TGeoHMatrix phiTransformation = currTransformation * TGeoTranslation(0, 0, zMin + (zMax - zMin) / 2.0);
    zoneList phiZone = cutByPhi(phiTransformation, sPhi, dPhi, rMaxMax, rMaxMax, (zMax - zMin) / 2.0);
//     startPoint = TGeoTranslation(currRotation * TGeoTranslation(0, 0, zMin + (zMax - zMin) / 2.0));
//     zoneList phiZone = cutByPhi(startPoint, currRotation, sPhi, dPhi, rMaxMax, rMaxMax, (zMax - zMin) / 2.0);
    out = andZone(out, phiZone);
    outerShell = andZone(outerShell, phiZone);
    return std::make_pair(out, outerShell);
}
inline std::pair<zoneList, zoneList> polyhedraToZones(TGeoPgon *polyhedra,TGeoHMatrix currTransformation, double scale) {
    zoneList out, outerShell; //outer shell for conjunction at generating mother volume
    TGeoHMatrix currRotation = TGeoRotation(currTransformation); //currRotation.SetTranslation(kNullVector);
    float sPhi = (float)polyhedra->GetPhi1(); //In Pcon angle converted to degrees, and TMatrixes work with degrees. (?)
    float dPhi = (float)polyhedra->GetDphi();
    int numSides = polyhedra->GetNedges();
    int nz = polyhedra->GetNz();
    float da = dPhi / numSides;
    float zCurrent, zPrev, rMaxCurrent, rMaxPrev, rMinCurrent, rMinPrev, rmin, rmax;
    double *zValues = polyhedra->GetZ();
    double *rMinValues = polyhedra->GetRmin();
    double *rMaxValues = polyhedra->GetRmax();
    float zMin = zValues[0] * scale, zMax = zValues[nz - 1] * scale;
    TGeoTranslation startZPoint, endZVec, vec0, vec1, vec2, vec3;
    TGeoHMatrix phiRotation;
    SGeoBody wed;
    zoneList currentZone;
    float zCenter = (zMax + zMin) / 2.0;
    //Work with external polyhedra
    for (int i = 1; i < nz; ++i) {
        zPrev = zCenter + (zValues[i - 1] - zCenter) * scale;
        zCurrent = zCenter + (zValues[i] - zCenter) * scale;
        if(zPrev == zCurrent) continue;
        rMaxPrev = rMaxValues[i - 1] * scale / cos(da/2*M_PI/180); //Convert angre to radians
        rMaxCurrent = rMaxValues[i] * scale / cos(da/2*M_PI/180);
        startZPoint = TGeoTranslation(currTransformation * TGeoTranslation(0, 0, zPrev));
        endZVec = TGeoTranslation(currRotation * TGeoTranslation(0, 0, zCurrent - zPrev));
        rmin = (rMaxCurrent >= rMaxPrev) ? rMaxPrev : rMaxCurrent;
        rmax = (rMaxCurrent >= rMaxPrev) ? rMaxCurrent : rMaxPrev;
        if (doubleEQ(rMaxPrev,0) && doubleEQ(rMaxCurrent,0)) continue;
        for (int k = 0; k < numSides; ++k) {
            phiRotation = TGeoHMatrix(); phiRotation.RotateZ(sPhi + da * k); phiRotation = currRotation * phiRotation;
            vec2 = TGeoTranslation(phiRotation * TGeoTranslation(1, 0, 0)); vec2 = vec2 * rmax;
            phiRotation = TGeoHMatrix(); phiRotation.RotateZ(sPhi + da * (k+1)); phiRotation = currRotation * phiRotation;
            vec3 = TGeoTranslation(phiRotation * TGeoTranslation(1, 0, 0)); vec3 = vec3 * rmax;
            wed = {1, {}};
            addVectorToElement(wed,0,startZPoint);
            addVectorToElement(wed,3,endZVec);
            addVectorToElement(wed,6,vec2);
            addVectorToElement(wed,9,vec3);
            currentZone = {{std::make_pair(1, wed)}};
            if (doubleNE(rMaxCurrent,rMaxPrev)) {
                vec0 = (rMaxCurrent < rMaxPrev) ? startZPoint + endZVec + vec3
                                                : startZPoint + vec2;
                vec1 = (rMaxCurrent < rMaxPrev) ? vec2 - vec3
                                                : vec3 - vec2;
                vec2 = (rMaxCurrent < rMaxPrev) ? endZVec * -1
                                                : endZVec;
                Double_t* a = (Double_t*)vec1.GetTranslation();
                Double_t* b = (Double_t*)vec2.GetTranslation();
                Double_t c[3] = { a[1] * b[2] - a[2] * b[1],
                                  a[2] * b[0] - a[0] * b[2],
                                  a[0] * b[1] - a[1] * b[0]
                                };
                Double_t vec3scale = (rmax-rmin) / sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);
                vec3 = TGeoTranslation(c[0], c[1], c[2]) * vec3scale;
                wed = {1, {}};
                addVectorToElement(wed,0,vec0);
                addVectorToElement(wed,3,vec1);
                addVectorToElement(wed,6,vec2);
                addVectorToElement(wed,9,vec3);
                currentZone = andZone(currentZone, {{std::make_pair(-1, wed)}});
            }
            outerShell = orZone(outerShell, currentZone);
        }
    }
    //Work with internal polyhedra
    for (int i = 1; i < nz; ++i) {
        zPrev = zCenter + (zValues[i - 1] - zCenter) * scale;
        zCurrent = zCenter + (zValues[i] - zCenter) * scale;
        if(zPrev == zCurrent) continue;
        rMinPrev = rMinValues[i - 1] * scale / cos(da/2*M_PI/180);
        rMinCurrent = rMinValues[i] * scale / cos(da/2*M_PI/180);
        startZPoint = TGeoTranslation(currRotation * TGeoTranslation(0, 0, zPrev));
        endZVec = TGeoTranslation(currRotation * TGeoTranslation(0, 0, zCurrent - zPrev));
        rmin = (rMinCurrent >= rMinPrev) ? rMinPrev : rMinCurrent;
        rmax = (rMinCurrent >= rMinPrev) ? rMinCurrent : rMinPrev;
        if (doubleEQ(rMinPrev,0) && doubleEQ(rMinCurrent,0)) continue;
        for (int k = 0; k < numSides; ++k) {
            phiRotation = TGeoHMatrix(); phiRotation.RotateZ(sPhi + da * k); phiRotation = currRotation * phiRotation;
            vec2 = TGeoTranslation(phiRotation * TGeoTranslation(1, 0, 0)); vec2 = vec2 * rmax;
            phiRotation = TGeoHMatrix(); phiRotation.RotateZ(sPhi + da * (k+1)); phiRotation = currRotation * phiRotation;
            vec3 = TGeoTranslation(phiRotation * TGeoTranslation(1, 0, 0)); vec3 = vec3 * rmax;
            wed = {1, {}};
            addVectorToElement(wed,0,startZPoint);
            addVectorToElement(wed,3,endZVec);
            addVectorToElement(wed,6,vec2);
            addVectorToElement(wed,9,vec3);
            currentZone = {{std::make_pair(1, wed)}};
            if (doubleNE(rMinCurrent,rMinPrev)) {
                vec0 = (rMinCurrent < rMinPrev) ? startZPoint + endZVec + vec3
                                                : startZPoint + vec2;
                vec1 = (rMinCurrent < rMinPrev) ? vec2 - vec3
                                                : vec3 - vec2;
                vec2 = (rMinCurrent < rMinPrev) ? endZVec * -1
                                                : endZVec;
                Double_t* a = (Double_t*)vec1.GetTranslation();
                Double_t* b = (Double_t*)vec2.GetTranslation();
                Double_t c[3] = { a[1] * b[2] - a[2] * b[1],
                                  a[2] * b[0] - a[0] * b[2],
                                  a[0] * b[1] - a[1] * b[0]
                                };
                Double_t vec3scale = (rmax-rmin) / sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);
                vec3 = TGeoTranslation(c[0], c[1], c[2]) * vec3scale;
                wed = {1, {}};
                addVectorToElement(wed,0,vec0);
                addVectorToElement(wed,3,vec1);
                addVectorToElement(wed,6,vec2);
                addVectorToElement(wed,9,vec3);
                currentZone = andZone(currentZone, {{std::make_pair(-1, wed)}});
            }
            out = orZone(out, currentZone);
        }
    }
    out = andZone(outerShell, notZone(out));
    return std::make_pair(out, outerShell);
}
inline std::pair<zoneList, zoneList> ellipticalTubeToZones(TGeoEltu *tube,TGeoHMatrix currTransformation, double scale) {
    zoneList out, outerShell; //outer shell for conjunction at generating mother volume
    TGeoHMatrix currRotation = TGeoHMatrix(currTransformation); currRotation.SetTranslation(kNullVector);
    float dx = (float)tube->GetA() * scale;
    float dy = (float)tube->GetB() * scale;
    float dz = (float)tube->GetDz() * scale;
    TGeoTranslation startPoint = TGeoTranslation(currTransformation * TGeoTranslation(0,0,-dz));
    TGeoTranslation endVec = TGeoTranslation(currRotation * TGeoTranslation(0,0,2*dz));
    TGeoTranslation xAxis = TGeoTranslation(currRotation * TGeoTranslation(dx,0,0));
    TGeoTranslation yAxis = TGeoTranslation(currRotation * TGeoTranslation(0,dy,0));
    SGeoBody rec = {6,{}};
    addVectorToElement(rec,0,startPoint);
    addVectorToElement(rec,3,endVec);
    addVectorToElement(rec,6,xAxis);
    addVectorToElement(rec,9,yAxis);
    out = {{std::make_pair(1,rec)}};
    outerShell = out;
    return std::make_pair(out, outerShell);
}
inline std::pair<zoneList, zoneList> arb8ToZones(TGeoArb8 *gtrap,TGeoHMatrix currTransformation, double scale) {
    zoneList out, outerShell; //outer shell for conjunction at generating mother volume
    TGeoHMatrix currRotation = TGeoHMatrix(currTransformation); currRotation.SetTranslation(kNullVector);
    float dz = (float)gtrap->GetDz() * scale;
    Double_t vertexes [8][2];
    Double_t* tmpPointer = gtrap->GetVertices();
    memcpy(&vertexes[0][0], tmpPointer, 8*2*sizeof(Double_t));
    TGeoTranslation tmp;
    float cx, cy, cz;
    SGeoBody trap = {2,{}};
    for(int i = 0; i<8; i++){
        cx = vertexes[i][0] * scale;
        cy = vertexes[i][1] * scale;
        cz = (i<4)?-dz:dz;
        tmp = TGeoTranslation(currTransformation * TGeoTranslation(cx,cy,cz));
        addVectorToElement(trap,i*3,tmp);
    }
    out = {{std::make_pair(1,trap)}};
    outerShell = out;
    return std::make_pair(out, outerShell);
}
std::pair<zoneList, zoneList> compositeShapeToZones(TGeoCompositeShape *shape,TGeoHMatrix currTransformation, double scale) {
    zoneList out, outerShell; //outer shell for conjunction at generating mother volume
    TGeoBoolNode* boolNode = shape->GetBoolNode();
    TGeoMatrix* lm = boolNode->GetLeftMatrix();
    TGeoShape* ls = boolNode->GetLeftShape();
    TGeoHMatrix ln = (*lm) * currTransformation;
    std::pair<zoneList, zoneList> lo = getZoneFromShape(ls,ln,scale);
    TGeoMatrix* rm = boolNode->GetRightMatrix();
    TGeoShape* rs = boolNode->GetRightShape();
    TGeoHMatrix rn = (*rm) * currTransformation;
    std::pair<zoneList, zoneList> ro = getZoneFromShape(rs,rn,scale);
    switch(boolNode->GetBooleanOperator()){
        case TGeoBoolNode::kGeoUnion:
            return std::pair<zoneList,zoneList>(orZone(lo.first,ro.first),orZone(lo.second,ro.second));
            break;
        case TGeoBoolNode::kGeoIntersection:
            return std::pair<zoneList,zoneList>(andZone(lo.first,ro.first),orZone(lo.second,ro.second));
            break;
        case TGeoBoolNode::kGeoSubtraction:
            return std::pair<zoneList,zoneList>(andZone(lo.first,ro.first),notZone(orZone(lo.second,ro.second)));
            break;
    }
}

std::pair<zoneList, MediumData> addOuterVacuum(TGeoNode *world) {
    zoneList out = getZoneFromBody(world, TGeoIdentity(), 1.5 * getDefaultScale()).first.second;
    zoneList outerWorld = getZoneFromBody(world, TGeoIdentity()).first.second;    
    MediumData medium = {1000, 0, 0, {}};
    return std::make_pair(andZone(out, notZone(outerWorld)), medium);
}

MediumData getMedium(TGeoNode *node) {
    double g = 1;
    double cm3 = getDefaultScale() * getDefaultScale() * getDefaultScale();
    MediumData out = {0, 0, 0, {}};
    TGeoMaterial *medium = node->GetMedium()->GetMaterial();
//     G4Material *mat = body->GetLogicalVolume()->GetMaterial();
    out.nChemEl = medium->GetNelements();
    out.nType = (out.nChemEl == 1)?1:2;
    out.Rho = medium->GetDensity() / (g/cm3);
    if (out.Rho == 0){
        printf("Is medium a Vacuum?");
        out.nType = 1000;
        return out;
    }
    for(int i = 0; i<out.nChemEl; ++i){
        out.Elements[i].Nuclid = medium->GetElement(i)->Z();
        out.Elements[i].A = medium->GetElement(i)->N();
        out.Elements[i].Z = medium->GetElement(i)->Z();
        out.Elements[i].Density = (medium->IsMixture()) ? out.Rho * (((TGeoMixture*)medium)->GetWmixt())[i] : out.Rho;
        out.Elements[i].Conc = out.Elements[i].Density / medium->GetElement(i)->A() * TMath::Na()* 1E-27;
//         out.Elements[i].Conc = mat->GetVecNbOfAtomsPerVolume()[i] * cm3 * 1E-27; //Correct
//         out.Elements[i].ionEv = mat->GetElement(i)->GetIonisation()->GetMeanExcitationEnergy() / eV; //As I unserstand, it is, but it can fill automatically.
        for(int k = 0; k<i; ++k){
            if(out.nType != 3 && out.Elements[i].Z == out.Elements[k].Z)
                out.nType = 3;
        }
    }
    return out;
}

}