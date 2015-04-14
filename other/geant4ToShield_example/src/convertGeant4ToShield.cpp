#include "convertGeant4ToShield.h"
namespace geant4toshield {
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
SGeoBody createCone(G4ThreeVector startPoint, G4ThreeVector vectorToEnd, float rStart, float rEnd) {
    SGeoBody outBody;
    if (rStart == rEnd) {
        outBody = {5, { (float)startPoint.x(), (float)startPoint.y(), (float)startPoint.z(),
                        (float)vectorToEnd.x(), (float)vectorToEnd.y(), (float)vectorToEnd.z(), rStart
                      }
                  };
    } else {
        startPoint = (rStart > rEnd) ? startPoint : startPoint + vectorToEnd;
        vectorToEnd = (rStart > rEnd) ? vectorToEnd : vectorToEnd * (-1);
        float rMax = (rStart > rEnd) ? rStart : rEnd;
        float rMin = (rStart > rEnd) ? rEnd : rStart;
        outBody = {7, { (float)startPoint.x(), (float)startPoint.y(), (float)startPoint.z(),
                        (float)vectorToEnd.x(), (float)vectorToEnd.y(), (float)vectorToEnd.z(), rMax, rMin
                      }
                  };
    }
    return outBody;
}
zoneList cutByPhi(G4ThreeVector currTranslation, G4RotationMatrix currRotation,
                  float sPhi, float dPhi, float halfX, float halfY, float halfZ) {
    zoneList outList;
    float ePhi = sPhi + dPhi;
    if(dPhi>=2*M_PI) return {{}};
    G4RotationMatrix innerRot = (G4RotationMatrix().rotateZ(sPhi)) * currRotation;
//         G4RotationMatrix innerRot = G4RotationMatrix(G4ThreeVector(0,0,1),sPhi) * currRotation;
    G4ThreeVector startFirstBox = currTranslation + innerRot * G4ThreeVector(0, -halfY, -halfZ);
    G4ThreeVector vec11 = innerRot * G4ThreeVector(-halfX, 0, 0);
    G4ThreeVector vec12 = innerRot * G4ThreeVector(0, 2 * halfY, 0);
    G4ThreeVector vec13 = innerRot * G4ThreeVector(0, 0, 2 * halfZ);
    SGeoBody box1 = {3, {(float)startFirstBox.x(), (float)startFirstBox.y(), (float)startFirstBox.z(),
                         (float)vec11.x(), (float)vec11.y(), (float)vec11.z(),
                         (float)vec12.x(), (float)vec12.y(), (float)vec12.z(),
                         (float)vec13.x(), (float)vec13.y(), (float)vec13.z()
                        }
                    };
    G4RotationMatrix outerRot = (G4RotationMatrix().rotateZ(ePhi)) * currRotation;
    G4ThreeVector startSecondBox = currTranslation + outerRot * G4ThreeVector(0, -halfY, -halfZ);
    G4ThreeVector vec21 = outerRot * G4ThreeVector(+halfX, 0, 0);
    G4ThreeVector vec22 = outerRot * G4ThreeVector(0, 2 * halfY, 0);
    G4ThreeVector vec23 = outerRot * G4ThreeVector(0, 0, 2 * halfZ);
    SGeoBody box2 = {3, {(float)startSecondBox.x(), (float)startSecondBox.y(), (float)startSecondBox.z(),
                         (float)vec21.x(), (float)vec21.y(), (float)vec21.z(),
                         (float)vec22.x(), (float)vec22.y(), (float)vec22.z(),
                         (float)vec23.x(), (float)vec23.y(), (float)vec23.z()
                        }
                    };
    if (dPhi <= M_PI) {
        outList = {{std::make_pair(-1, box1), std::make_pair(-1, box2)}};
    } else {
        outList = {{std::make_pair(-1, box1)}, {std::make_pair(-1, box2)}};
    }
    return outList;
}

inline std::pair<zoneList, zoneList> boxToZones(G4Box *box, G4ThreeVector currTranslation, G4RotationMatrix currRotation, double scale) {
    zoneList out, outerShell; //outer shell for conjunction at generating mother volume
    float x = (float) box->GetXHalfLength() * scale;
    float y = (float) box->GetYHalfLength() * scale;
    float z = (float) box->GetZHalfLength() * scale;
    G4ThreeVector startVector = currTranslation + currRotation * G4ThreeVector(-x, -y, -z);
    G4ThreeVector vec1 = currRotation * G4ThreeVector(2 * x, 0, 0);
    G4ThreeVector vec2 = currRotation * G4ThreeVector(0, 2 * y, 0);
    G4ThreeVector vec3 = currRotation * G4ThreeVector(0, 0, 2 * z);
    SGeoBody tmp = {3, { (float)startVector.x(), (float)startVector.y(), (float)startVector.z(),
                         (float)vec1.x(), (float)vec1.y(), (float)vec1.z(),
                         (float)vec2.x(), (float)vec2.y(), (float)vec2.z(),
                         (float)vec3.x(), (float)vec3.y(), (float)vec3.z()
                       }
                   };
//         out.push_back(std::vector<zoneElement>(1,std::make_pair(1,tmp)));
    out = {{std::make_pair(1, tmp)}};
//         outerShell.push_back(std::vector<zoneElement>(1,std::make_pair(1,tmp)));
    outerShell = {{std::make_pair(1, tmp)}};
    return std::make_pair(out, outerShell);
}
inline std::pair<zoneList, zoneList> tubeToZones(G4Tubs *tube, G4ThreeVector currTranslation, G4RotationMatrix currRotation, double scale) {
    zoneList out, outerShell; //outer shell for conjunction at generating mother volume
    float rIn = (float) tube->GetInnerRadius() * scale;
    float rOut = (float) tube->GetOuterRadius() * scale;
    float z = (float) tube->GetZHalfLength() * scale;
    float sPhi = (float) tube->GetStartPhiAngle(); // or tube->GetSPhi() ?
    float dPhi = (float) tube->GetDeltaPhiAngle(); // or tube->GetDPhi() ?
    G4ThreeVector startTube = currTranslation + currRotation * G4ThreeVector(0, 0, -z);
    G4ThreeVector endTubeVec = currRotation * G4ThreeVector(0, 0, 2 * z);
    SGeoBody innerTube = {5, {(float)startTube.x(), (float)startTube.y(), (float)startTube.z(),
                              (float)endTubeVec.x(), (float)endTubeVec.y(), (float)endTubeVec.z(), rIn
                             }
                         };
    SGeoBody outerTube = {5, {(float)startTube.x(), (float)startTube.y(), (float)startTube.z(),
                              (float)endTubeVec.x(), (float)endTubeVec.y(), (float)endTubeVec.z(), rOut
                             }
                         };

    out = {{std::make_pair(1, outerTube)}};
    outerShell = {{std::make_pair(1, outerTube)}};
    if (rIn != 0) {
        out = andZone(out, {{std::make_pair(-1, innerTube)}});
    }
    zoneList phiZone = cutByPhi(currTranslation, currRotation, sPhi, dPhi, rOut, rOut, z);
    out = andZone(out, phiZone);
    outerShell = andZone(outerShell, phiZone);
    return std::make_pair(out, outerShell);
}
inline std::pair<zoneList, zoneList> coneToZones(G4Cons *cone, G4ThreeVector currTranslation, G4RotationMatrix currRotation, double scale) {
    zoneList out, outerShell; //outer shell for conjunction at generating mother volume
    float rMin1 = (float)cone->GetInnerRadiusMinusZ() * scale;
    float rMin2 = (float)cone->GetInnerRadiusPlusZ() * scale;
    float rMax1 = (float)cone->GetOuterRadiusMinusZ() * scale;
    float rMax2 = (float)cone->GetOuterRadiusPlusZ() * scale;
    float z = (float)cone->GetZHalfLength() * scale;
    float sPhi = (float) cone->GetStartPhiAngle();
    float dPhi = (float) cone->GetDeltaPhiAngle();
    G4ThreeVector startTubeVector = currTranslation + currRotation * G4ThreeVector(0, 0, -z);
    G4ThreeVector endTubeVec = currRotation * G4ThreeVector(0, 0, 2 * z);
    SGeoBody outerCone = createCone(startTubeVector, endTubeVec, rMax1, rMax2);
    SGeoBody innerCone = createCone(startTubeVector, endTubeVec, rMin1, rMin2);

    out = {{std::make_pair(1, outerCone)}};
    outerShell = {{std::make_pair(1, outerCone)}};
    if (rMin1 != 0 || rMin2 != 0) {
        out = andZone(out, {{std::make_pair(-1, innerCone)}});
    }

    float rOut = (rMax1 > rMax2) ? rMax1 : rMax2;
    zoneList phiZone = cutByPhi(currTranslation, currRotation, sPhi, dPhi, rOut, rOut, z);
    out = andZone(out, phiZone);
    outerShell = andZone(outerShell, phiZone);

    return std::make_pair(out, outerShell);
}
inline std::pair<zoneList, zoneList> trdToZones(G4Trd *trd, G4ThreeVector currTranslation, G4RotationMatrix currRotation, double scale) {
    zoneList out, outerShell; //outer shell for conjunction at generating mother volume
    float dz = (float) trd->GetZHalfLength() * scale;
    float dx1 = (float) trd->GetXHalfLength1() * scale;
    float dy1 = (float) trd->GetYHalfLength1() * scale;
    float dx2 = (float) trd->GetXHalfLength2() * scale;
    float dy2 = (float) trd->GetYHalfLength2() * scale;
    G4ThreeVector v1 = currTranslation + currRotation * G4ThreeVector(-dx1, -dy1, -dz);
    G4ThreeVector v2 = currTranslation + currRotation * G4ThreeVector(-dx1, dy1, -dz);
    G4ThreeVector v3 = currTranslation + currRotation * G4ThreeVector(dx1, dy1, -dz);
    G4ThreeVector v4 = currTranslation + currRotation * G4ThreeVector(dx1, -dy1, -dz);
    G4ThreeVector v5 = currTranslation + currRotation * G4ThreeVector(-dx2, -dy2, dz);
    G4ThreeVector v6 = currTranslation + currRotation * G4ThreeVector(-dx2, dy2, dz);
    G4ThreeVector v7 = currTranslation + currRotation * G4ThreeVector(dx2, dy2, dz);
    G4ThreeVector v8 = currTranslation + currRotation * G4ThreeVector(dx2, -dy2, dz);
    SGeoBody tmp = {2, { (float)v1.x(), (float)v1.y(), (float)v1.z(),
                         (float)v2.x(), (float)v2.y(), (float)v2.z(),
                         (float)v3.x(), (float)v3.y(), (float)v3.z(),
                         (float)v4.x(), (float)v4.y(), (float)v4.z(),
                         (float)v5.x(), (float)v5.y(), (float)v5.z(),
                         (float)v6.x(), (float)v6.y(), (float)v6.z(),
                         (float)v7.x(), (float)v7.y(), (float)v7.z(),
                         (float)v8.x(), (float)v8.y(), (float)v8.z()
                       }
                   };
    out = {{std::make_pair(1, tmp)}};
    outerShell = out;
    return std::make_pair(out, outerShell);
}
inline std::pair<zoneList, zoneList> sphereToZones(G4Sphere *sphere, G4ThreeVector currTranslation, G4RotationMatrix currRotation, double scale) {
    zoneList out, outerShell; //outer shell for conjunction at generating mother volume

    float rIn = (float) sphere->GetInnerRadius() * scale;
    float rOut = (float) sphere->GetOuterRadius() * scale;
    float sPhi = (float) sphere->GetStartPhiAngle();
    float dPhi = (float) sphere->GetDeltaPhiAngle();
    float sTheta = (float) sphere->GetStartThetaAngle();
    float dTheta = (float) sphere->GetDeltaThetaAngle();
    float eTheta = sTheta + dTheta;
    SGeoBody innerSphere = {0, {(float)currTranslation.x(), (float)currTranslation.y(), (float)currTranslation.z(), rIn}};
    SGeoBody outerSphere = {0, {(float)currTranslation.x(), (float)currTranslation.y(), (float)currTranslation.z(), rOut}};

    out = {{std::make_pair(1, outerSphere)}};
    outerShell = {{std::make_pair(1, outerSphere)}};
    if (rIn != 0) {
        out = andZone(out, {{std::make_pair(-1, innerSphere)}});
    }
    zoneList phiZone = cutByPhi(currTranslation, currRotation, sPhi, dPhi, rOut, rOut, rOut);
    out = andZone(out, phiZone);
    outerShell = andZone(outerShell, phiZone);

    if (sTheta != 0 && sTheta != M_PI && sTheta != (M_PI / 2.0)) {
        G4ThreeVector trc1Start = currTranslation + currRotation * G4ThreeVector(0, 0, -rOut * cos(sTheta));
        G4ThreeVector trc1Vec = currRotation * G4ThreeVector(0, 0, (rOut - rIn) * cos(sTheta));
        SGeoBody trc1 = {7, {(float)trc1Start.x(), (float)trc1Start.y(), (float)trc1Start.z(),
                             (float)trc1Vec.x(), (float)trc1Vec.y(), (float)trc1Vec.z(),
                             (float)(rOut * tan(sTheta)), (float)(rIn * sin(sTheta))
                            }
                        };
        int mul = (sTheta < M_PI / 2.0) ? -1 : 1;
        out = andZone(out, {{std::make_pair(mul, trc1)}});
        outerShell = andZone(outerShell, {{std::make_pair(mul, trc1)}});
    } else if (sTheta == (M_PI / 2.0)) {
        G4ThreeVector db1s = currTranslation + currRotation * G4ThreeVector(-rOut, -rOut, 0);
        G4ThreeVector db1v1 = currRotation * G4ThreeVector(2 * rOut, 0, 0);
        G4ThreeVector db1v2 = currRotation * G4ThreeVector(0, 2 * rOut, 0);
        G4ThreeVector db1v3 = currRotation * G4ThreeVector(0, 0, -rOut);
        SGeoBody db1 = {3, {(float)db1s.x(), (float)db1s.y(), (float)db1s.z(),
                            (float)db1v1.x(), (float)db1v1.y(), (float)db1v1.z(),
                            (float)db1v2.x(), (float)db1v2.y(), (float)db1v2.z(),
                            (float)db1v3.x(), (float)db1v3.y(), (float)db1v3.z()
                           }
                       };
        out = andZone(out, {{std::make_pair(-1, db1)}});
        outerShell = andZone(outerShell, {{std::make_pair(-1, db1)}});
    }

    if (eTheta != 0 && eTheta != M_PI && eTheta != (M_PI / 2.0)) {
        G4ThreeVector trc2Start = currTranslation + currRotation * G4ThreeVector(0, 0, -rOut * cos(eTheta));
        G4ThreeVector trc2Vec = currRotation * G4ThreeVector(0, 0, (rOut - rIn) * cos(eTheta));
        SGeoBody trc2 = {7, {(float)trc2Start.x(), (float)trc2Start.y(), (float)trc2Start.z(),
                             (float)trc2Vec.x(), (float)trc2Vec.y(), (float)trc2Vec.z(),
                             (float)(rOut * sin(eTheta)), (float)(rIn * tan(eTheta))
                            }
                        };
        int mul = (eTheta > M_PI / 2.0) ? -1 : 1;
        out = andZone(out, {{std::make_pair(mul, trc2)}});
        outerShell = andZone(outerShell, {{std::make_pair(mul, trc2)}});
    } else if (eTheta == (M_PI / 2.0)) {
        G4ThreeVector db2s = currTranslation + currRotation * G4ThreeVector(-rOut, -rOut, 0);
        G4ThreeVector db2v1 = currRotation * G4ThreeVector(2 * rOut, 0, 0);
        G4ThreeVector db2v2 = currRotation * G4ThreeVector(0, 2 * rOut, 0);
        G4ThreeVector db2v3 = currRotation * G4ThreeVector(0, 0, rOut);
        SGeoBody db2 = {3, {(float)db2s.x(), (float)db2s.y(), (float)db2s.z(),
                            (float)db2v1.x(), (float)db2v1.y(), (float)db2v1.z(),
                            (float)db2v2.x(), (float)db2v2.y(), (float)db2v2.z(),
                            (float)db2v3.x(), (float)db2v3.y(), (float)db2v3.z()
                           }
                       };
        out = andZone(out, {{std::make_pair(-1, db2)}});
        outerShell = andZone(outerShell, {{std::make_pair(-1, db2)}});
    }
    return std::make_pair(out, outerShell);
}
inline std::pair<zoneList, zoneList> orbToZones(G4Orb *orb, G4ThreeVector currTranslation, G4RotationMatrix currRotation, double scale) {
    zoneList out, outerShell; //outer shell for conjunction at generating mother volume
    float r = (float) orb->GetRadius() * scale;
    SGeoBody sph = {0, {(float)currTranslation.x(), (float)currTranslation.y(), (float)currTranslation.z(), r}};
    out = {{std::make_pair(1, sph)}};
    outerShell = out;
    return std::make_pair(out, outerShell);
}
inline std::pair<zoneList, zoneList> polyconeToZones(G4Polycone *polycone, G4ThreeVector currTranslation, G4RotationMatrix currRotation, double scale) {
    zoneList out, outerShell; //outer shell for conjunction at generating mother volume
    G4PolyconeHistorical *hist = polycone->GetOriginalParameters();
    float sPhi = (float)hist->Start_angle;
    float dPhi = (float)hist->Opening_angle;
    int numPlanes = hist->Num_z_planes;
    float zCurrent, zPrev, rMaxCurrent, rMaxPrev, rMinCurrent, rMinPrev;
    float zMin = hist->Z_values[0], zMax = hist->Z_values[numPlanes - 1], rMaxMax = hist->Rmax[0];
    G4ThreeVector startPoint, endVec;
    zoneElement currentOuterCone, currentInnerCone;
    float zCenter = (zMax + zMin) / 2.0;
    for (int i = 1; i < numPlanes; ++i) {
        zPrev = zCenter + (hist->Z_values[i - 1] - zCenter) * scale;
        zCurrent = zCenter + (hist->Z_values[i] - zCenter) * scale;
        rMaxPrev = hist->Rmax[i - 1] * scale;
        rMaxCurrent = hist->Rmax[i] * scale;
        rMinPrev = hist->Rmin[i - 1] * scale;
        rMinCurrent = hist->Rmin[i] * scale;
        startPoint = currTranslation + currRotation * G4ThreeVector(0, 0, zPrev);
        endVec = currRotation * G4ThreeVector(0, 0, zCurrent - zPrev);
        rMaxMax = (rMaxCurrent > rMaxMax) ? rMaxCurrent : rMaxMax;
        if (rMaxPrev != 0 || rMaxCurrent != 0) {
            currentOuterCone = std::make_pair(1, createCone(startPoint, endVec, rMaxPrev, rMaxCurrent));
            outerShell = orZone(outerShell, {{currentOuterCone}});
        }
        if (rMinPrev != 0 || rMinCurrent != 0) {
            currentInnerCone = std::make_pair(1, createCone(startPoint, endVec, rMinPrev, rMinCurrent));
            out = orZone(out, {{currentInnerCone}});
        }
    }
    out = andZone(outerShell, notZone(out));
    zoneList phiZone = cutByPhi(currTranslation + currRotation * G4ThreeVector(0, 0, zMin + (zMax - zMin) / 2.0),
                                currRotation, sPhi, dPhi, rMaxMax, rMaxMax, (zMax - zMin) / 2.0);
    out = andZone(out, phiZone);
    outerShell = andZone(outerShell, phiZone);
    return std::make_pair(out, outerShell);
}
inline std::pair<zoneList, zoneList> polyhedraToZones(G4Polyhedra *polyhedra, G4ThreeVector currTranslation, G4RotationMatrix currRotation, double scale) {
    zoneList out, outerShell; //outer shell for conjunction at generating mother volume
    G4PolyhedraHistorical *hist = polyhedra->GetOriginalParameters();
    float sPhi = (float)hist->Start_angle;
    float dPhi = (float)hist->Opening_angle;
    int numSides = hist->numSide;
    int numPlanes = hist->Num_z_planes;
    float da = dPhi / numSides;
    float zCurrent, zPrev, rMaxCurrent, rMaxPrev, rMinCurrent, rMinPrev, rmin, rmax;
    float zMin = hist->Z_values[0], zMax = hist->Z_values[numPlanes - 1], rMaxMax = hist->Rmax[0] / cos(da / 2.0);
    G4ThreeVector startZPoint, endZVec, vec0, vec1, vec2, vec3;
    SGeoBody wed;
    zoneList currentZone;
    float zCenter = (zMax + zMin) / 2.0;
    //Work with external polycone
    for (int i = 1; i < numPlanes; ++i) {
        zPrev = zCenter + (hist->Z_values[i - 1] - zCenter) * scale;
        zCurrent = zCenter + (hist->Z_values[i] - zCenter) * scale;
        rMaxPrev = hist->Rmax[i - 1] * scale;
        rMaxCurrent = hist->Rmax[i] * scale;
        startZPoint = currTranslation + currRotation * G4ThreeVector(0, 0, zPrev);
        endZVec = currRotation * G4ThreeVector(0, 0, zCurrent - zPrev);
        rMaxMax = (rMaxCurrent / cos(da / 2.0) > rMaxMax) ? rMaxCurrent / cos(da / 2.0) : rMaxMax;
        rmin = (rMaxCurrent >= rMaxPrev) ? rMaxCurrent : rMaxPrev;
        rmax = (rMaxCurrent >= rMaxPrev) ? rMaxCurrent : rMaxPrev;
        if (rMaxPrev != 0 || rMaxCurrent != 0) {
            for (int k = 0; i < numSides; ++k) {
                vec2 = currRotation * G4ThreeVector(1, 0, 0).rotateZ(sPhi + da * k) * rmax;
                vec3 = currRotation * G4ThreeVector(1, 0, 0).rotateZ(sPhi + da * (k + 1)) * rmax;
                wed = {1, { (float)startZPoint.x(), (float)startZPoint.y(), (float)startZPoint.z(),
                            (float)endZVec.x(), (float)endZVec.y(), (float)endZVec.z(),
                            (float)vec2.x(), (float)vec2.y(), (float)vec2.z(),
                            (float)vec3.x(), (float)vec3.y(), (float)vec3.z()
                          }
                      };
                currentZone = {{std::make_pair(1, wed)}};
                if (rMaxCurrent != rMaxPrev) {
                    vec0 = startZPoint + vec2;
                    vec1 = vec3 - vec2;
                    vec2 = (rMaxCurrent < rMaxPrev) ? endZVec + vec2 * (rmin / rmax) - vec2
                           : endZVec + vec2 - vec2 * (rmin / rmax);
                    vec3 = (rMaxCurrent < rMaxPrev) ? endZVec
                           : -endZVec;
                    wed = {1, { (float)vec0.x(), (float)vec0.y(), (float)vec0.z(),
                                (float)vec1.x(), (float)vec1.y(), (float)vec1.z(),
                                (float)vec2.x(), (float)vec2.y(), (float)vec2.z(),
                                (float)vec3.x(), (float)vec3.y(), (float)vec3.z()
                              }
                          };
                    currentZone = andZone(currentZone, {{std::make_pair(-1, wed)}});
                }
                outerShell = orZone(outerShell, currentZone);
            }
        }
    }
    //Work with internal polycone
    for (int i = 1; i < numPlanes; ++i) {
        zPrev = zCenter + (hist->Z_values[i - 1] - zCenter) * scale;
        zCurrent = zCenter + (hist->Z_values[i] - zCenter) * scale;
        rMinPrev = hist->Rmin[i - 1] * scale;
        rMinCurrent = hist->Rmin[i] * scale;
        startZPoint = currTranslation + currRotation * G4ThreeVector(0, 0, zPrev);
        endZVec = currRotation * G4ThreeVector(0, 0, zCurrent - zPrev);
        rmin = (rMinCurrent >= rMinPrev) ? rMinCurrent : rMinPrev;
        rmax = (rMinCurrent >= rMinPrev) ? rMinCurrent : rMinPrev;
        if (rMinPrev != 0 || rMinCurrent != 0) {
            for (int k = 0; i < numSides; ++k) {
                vec2 = currRotation * G4ThreeVector(1, 0, 0).rotateZ(sPhi + da * k) * rmax;
                vec3 = currRotation * G4ThreeVector(1, 0, 0).rotateZ(sPhi + da * (k + 1)) * rmax;
                wed = {1, { (float)startZPoint.x(), (float)startZPoint.y(), (float)startZPoint.z(),
                            (float)endZVec.x(), (float)endZVec.y(), (float)endZVec.z(),
                            (float)vec2.x(), (float)vec2.y(), (float)vec2.z(),
                            (float)vec3.x(), (float)vec3.y(), (float)vec3.z()
                          }
                      };
                currentZone = {{std::make_pair(1, wed)}};
                if (rMinCurrent != rMinPrev) {
                    vec0 = startZPoint + vec2;
                    vec1 = vec3 - vec2;
                    vec2 = (rMinCurrent < rMinPrev) ? endZVec + vec2 * (rmin / rmax) - vec2
                           : endZVec + vec2 - vec2 * (rmin / rmax);
                    vec3 = (rMinCurrent < rMinPrev) ? endZVec
                           : -endZVec;
                    wed = {1, { (float)vec0.x(), (float)vec0.y(), (float)vec0.z(),
                                (float)vec1.x(), (float)vec1.y(), (float)vec1.z(),
                                (float)vec2.x(), (float)vec2.y(), (float)vec2.z(),
                                (float)vec3.x(), (float)vec3.y(), (float)vec3.z()
                              }
                          };
                    currentZone = andZone(currentZone, {{std::make_pair(-1, wed)}});
                }
                out = orZone(out, currentZone);
            }
        }
    }
    out = andZone(outerShell, notZone(out));
    return std::make_pair(out, outerShell);
}
inline std::pair<zoneList, zoneList> ellipticalTubeToZones(G4EllipticalTube *tube, G4ThreeVector currTranslation, G4RotationMatrix currRotation, double scale) {
    zoneList out, outerShell; //outer shell for conjunction at generating mother volume
    float dx = (float)tube->GetDx() * scale;
    float dy = (float)tube->GetDy() * scale;
    float dz = (float)tube->GetDz() * scale;
    G4ThreeVector startPoint = currTranslation + currRotation * G4ThreeVector(0,0,-dz);
    G4ThreeVector endVec = currRotation * G4ThreeVector(0,0,2*dz);
    G4ThreeVector xAxis = currRotation * G4ThreeVector(dx,0,0);
    G4ThreeVector yAxis = currRotation * G4ThreeVector(0,dy,0);
    SGeoBody rec = {6,{(float)startPoint.x(), (float)startPoint.y(), (float)startPoint.z(),
                        (float)endVec.x(), (float)endVec.y(), (float)endVec.z(),
                        (float)xAxis.x(), (float)xAxis.y(), (float)xAxis.z(),
                        (float)yAxis.x(), (float)yAxis.y(), (float)yAxis.z() }};
    out = {{std::make_pair(1,rec)}};
    outerShell = out;
    return std::make_pair(out, outerShell);
}
inline std::pair<zoneList, zoneList> ellipsoidToZones(G4Ellipsoid *ell, G4ThreeVector currTranslation, G4RotationMatrix currRotation, double scale) {
    zoneList out, outerShell; //outer shell for conjunction at generating mother volume
    float dx = (float) ell->GetSemiAxisMax(0) * scale;
    float dy = (float) ell->GetSemiAxisMax(1) * scale;
    float dz = (float) ell->GetSemiAxisMax(2) * scale;
    float bCut = (float) ell->GetZBottomCut() * scale;
    float tCut = (float) ell->GetZTopCut() * scale;
    
    G4ThreeVector xAxis = currRotation * G4ThreeVector(dx,0,0);
    G4ThreeVector yAxis = currRotation * G4ThreeVector(0,dy,0);
    G4ThreeVector zAxis = currRotation * G4ThreeVector(0,0,dz);
    SGeoBody ellBody = {8,{(float)currTranslation.x(), (float)currTranslation.y(), (float)currTranslation.z(),
                        (float)xAxis.x(), (float)xAxis.y(), (float)xAxis.z(),
                        (float)yAxis.x(), (float)yAxis.y(), (float)yAxis.z(),
                        (float)zAxis.x(), (float)zAxis.y(), (float)zAxis.z() }};
    out = {{std::make_pair(1,ellBody)}};
    G4ThreeVector xBox = currRotation * G4ThreeVector(2 * dx,0,0);
    G4ThreeVector yBox = currRotation * G4ThreeVector(0, 2 * dy,0);
    G4ThreeVector startBox, zBox ;// = currRotation * G4ThreeVector(0, 0,2 * dz);
    SGeoBody cutBox;
    if(bCut > (-dz)){
        startBox = currTranslation + currRotation * G4ThreeVector(-dx,-dy,bCut);
        zBox = currRotation * G4ThreeVector(0, 0, -dz - bCut);
        cutBox = {3, { (float)startBox.x(),(float)startBox.y(),(float)startBox.z(),
                        (float)xBox.x(), (float)xBox.y(), (float)xBox.z(),
                        (float)yBox.x(), (float)yBox.y(), (float)yBox.z(),
                        (float)zBox.x(), (float)zBox.y(), (float)zBox.z() }};
        out = andZone(out,{{std::make_pair(-1,cutBox)}});
    }
    if(tCut < dz){
        startBox = currTranslation + currRotation * G4ThreeVector(-dx,-dy,tCut);
        zBox = currRotation * G4ThreeVector(0, 0, dz - tCut);
        cutBox = {3, { (float)startBox.x(),(float)startBox.y(),(float)startBox.z(),
                        (float)xBox.x(), (float)xBox.y(), (float)xBox.z(),
                        (float)yBox.x(), (float)yBox.y(), (float)yBox.z(),
                        (float)zBox.x(), (float)zBox.y(), (float)zBox.z() }};
        out = andZone(out,{{std::make_pair(-1,cutBox)}});
    }
    outerShell = out;
    return std::make_pair(out, outerShell);
}
inline std::pair<zoneList, zoneList> genericTrapTubeToZones(G4GenericTrap *gtrap, G4ThreeVector currTranslation, G4RotationMatrix currRotation, double scale) {
    zoneList out, outerShell; //outer shell for conjunction at generating mother volume
    float dz = (float)gtrap->GetZHalfLength() * scale;
    std::vector<G4TwoVector> vertexes = gtrap->GetVertices();
    G4ThreeVector tmp;
    float cx, cy, cz;
    SGeoBody trap = {2,{}};
    for(int i = 0; i<8; i++){
        cx = vertexes.at(i).x() * scale;
        cy = vertexes.at(i).y() * scale;
        cz = (i<4)?-dz:dz;
        tmp = currTranslation + currRotation * G4ThreeVector(cx,cy,cz);
        trap.parameters[i*3] = tmp.x();
        trap.parameters[i*3+1] = tmp.y();
        trap.parameters[i*3+2] = tmp.z();        
    }
    out = {{std::make_pair(1,trap)}};
    outerShell = out;
    return std::make_pair(out, outerShell);
}

std::pair<std::pair<zoneList, zoneList>, std::pair<G4RotationMatrix, G4ThreeVector> >
getZoneFromBody(G4VPhysicalVolume *physBody, G4RotationMatrix oldRotation, G4ThreeVector oldTranslation, double scale) {//TODO
//  WARNING: extended initializer lists only available with -std=c++11 or -std=gnu++11
    std::pair<zoneList, zoneList> zonePair;
//     zoneList out, outerShell; //outer shell for conjunction at generating mother volume
    G4ThreeVector newTranslation = physBody->GetObjectTranslation() * scale;
    G4ThreeVector currTranslation = oldTranslation +  oldRotation * newTranslation;
    G4RotationMatrix newRotation = physBody->GetObjectRotationValue();
    G4RotationMatrix currRotation = newRotation * oldRotation;
    G4LogicalVolume *logicalBody = physBody->GetLogicalVolume();
    G4VSolid *solidBody = logicalBody->GetSolid();
    if (dynamic_cast<G4Box *>(solidBody) != NULL) {
        zonePair = boxToZones((G4Box *)solidBody, currTranslation, currRotation, scale);
    } else if (dynamic_cast<G4Tubs *>(solidBody) != NULL) {
        zonePair = tubeToZones((G4Tubs *)solidBody, currTranslation, currRotation, scale);
    } else if (dynamic_cast<G4CutTubs *>(solidBody) != NULL) {
        printf("G4CutTubs not implemented yet\n");
    } else if (dynamic_cast<G4Cons *>(solidBody) != NULL) {
        zonePair = coneToZones((G4Cons *)solidBody, currTranslation, currRotation, scale);
    } else if (dynamic_cast<G4Para *>(solidBody) != NULL) {
        printf("G4Para not implemented yet\n");
    } else if (dynamic_cast<G4Trd *>(solidBody) != NULL) {
        zonePair = trdToZones((G4Trd *)solidBody, currTranslation, currRotation, scale);
    } else if (dynamic_cast<G4Trap *>(solidBody) != NULL) { //I connot understand angles and parameters from getters
        printf("G4Trap not implemented yet\n");
    } else if (dynamic_cast<G4Orb *>(solidBody) != NULL) {
        zonePair = orbToZones((G4Orb *)solidBody, currTranslation, currRotation, scale);
    } else if (dynamic_cast<G4Sphere *>(solidBody) != NULL) {
        zonePair = sphereToZones((G4Sphere *)solidBody, currTranslation, currRotation, scale);
    } else if (dynamic_cast<G4Polycone *>(solidBody) != NULL) {
        zonePair = polyconeToZones((G4Polycone *)solidBody, currTranslation, currRotation, scale);
    } else if (dynamic_cast<G4Polyhedra *>(solidBody) != NULL) {
        zonePair = polyhedraToZones((G4Polyhedra *)solidBody, currTranslation, currRotation, scale);
    } else if (dynamic_cast<G4EllipticalTube *>(solidBody) != NULL) {
        zonePair = ellipticalTubeToZones((G4EllipticalTube *)solidBody, currTranslation, currRotation, scale);
    } else if (dynamic_cast<G4Ellipsoid*>(solidBody) != NULL) {
        zonePair = ellipsoidToZones((G4Ellipsoid*)solidBody, currTranslation, currRotation, scale);
    } else if (dynamic_cast<G4GenericTrap*>(solidBody) != NULL) {
        zonePair = genericTrapTubeToZones((G4GenericTrap*)solidBody, currTranslation, currRotation, scale);
    } else if (dynamic_cast<G4Tet*>(solidBody) != NULL) {
        printf("G4Tet not implemented yet\n");
    } else {
        printf("ELSE");
    }
    return std::make_pair(zonePair, std::make_pair(currRotation, currTranslation));
}
std::pair<zoneList, MediumData> addOuterVacuum(G4VPhysicalVolume *world) {
    zoneList out = getZoneFromBody(world, G4RotationMatrix(), G4ThreeVector(), 1.5 * getDefaultScale()).first.second;
    zoneList outerWorld = getZoneFromBody(world, G4RotationMatrix(), G4ThreeVector()).first.first;    
    MediumData medium = {0, 0, 0, {}};
    return std::make_pair(andZone(out, notZone(outerWorld)), medium);
}
MediumData getMedium(G4VPhysicalVolume *body) {
    MediumData out = {0, 0, 0, {}};
    G4Material *mat = body->GetLogicalVolume()->GetMaterial();
    out.nChemEl = mat->GetNumberOfElements();
    out.nType = (out.nChemEl == 1)?1:2;
    out.Rho = mat->GetDensity() / (g/cm3);
    for(int i = 0; i<out.nChemEl; ++i){
        G4IonisParamElm ip =*(mat->GetElement(i)->GetIonisation());
        out.Elements[i].Nuclid = mat->GetElement(i)->GetZ(); //Correct ?
        out.Elements[i].A = mat->GetElement(i)->GetN(); //Correct
        out.Elements[i].Z = mat->GetElement(i)->GetZ(); //Correct
        out.Elements[i].Density = out.Rho * mat->GetFractionVector()[i]; //Correct ?
        out.Elements[i].Conc = mat->GetVecNbOfAtomsPerVolume()[i] * cm3 * 1E-27; //Correct
        out.Elements[i].ionEv = mat->GetElement(i)->GetIonisation()->GetMeanExcitationEnergy() / eV; //As I unserstand, it is, but it can fill automatically.
        for(int k = 0; k<i; ++k){
            if(out.nType != 3 && out.Elements[i].A == out.Elements[k].A)
                out.nType = 3;
        }
    }
    return out;
}

std::pair < zoneList, std::pair<G4RotationMatrix, G4ThreeVector> >
convertToSGeoZone(G4VPhysicalVolume *physBody, G4RotationMatrix oldRotation, G4ThreeVector oldTranslation) { //TODO
    std::pair<std::pair<zoneList, zoneList>, std::pair<G4RotationMatrix, G4ThreeVector> > tmp = getZoneFromBody(physBody, oldRotation, oldTranslation);
    zoneList outZone = tmp.first.first;
    G4LogicalVolume *logicalBody = physBody->GetLogicalVolume();
    int countOfDaughters = logicalBody->GetNoDaughters();
    G4VPhysicalVolume *daughter;
    for (int i = 0; i < countOfDaughters; ++i) {
        daughter = logicalBody->GetDaughter(i);
        std::pair<zoneList, zoneList> daughterZone = getZoneFromBody(daughter, tmp.second.first, tmp.second.second).first;
        outZone = andZone(outZone, notZone(daughterZone.second)); //Only outel shell of body.
    }
    return std::make_pair(outZone, tmp.second);
}

zoneData convertBodyRecursive(G4VPhysicalVolume *physBody, G4RotationMatrix oldRotation, G4ThreeVector oldTranslation) {
    zoneData out;
    std::pair<zoneList, std::pair<G4RotationMatrix, G4ThreeVector> > tmp = convertToSGeoZone(physBody, oldRotation, oldTranslation);
    out.push_back(std::make_pair(tmp.first, getMedium(physBody)));
    G4LogicalVolume *logicalBody = physBody->GetLogicalVolume();
    int countOfDaughters = logicalBody->GetNoDaughters();
    G4VPhysicalVolume *daughter;
    zoneData daughterZones;
    for (int i = 0; i < countOfDaughters; ++i) {
        daughter = logicalBody->GetDaughter(i);
        zoneData daughterZones = convertBodyRecursive(daughter, tmp.second.first, tmp.second.second);
        out.insert(out.end(), daughterZones.begin(), daughterZones.end());
    }
    return out;
}
Geant4ShieldData zoneDataToShield(zoneData data) {
    Geant4ShieldData out;
    Geant4ShieldElement tmpElement;
    std::pair<zoneList, MediumData> currentPair;
    zoneList currentList;
    std::vector<SGeoBody> bodyVector;
    std::vector<int> zoneVector;
    for (unsigned int kk = 0; kk < data.size(); ++kk) {
        bodyVector.clear();
        zoneVector.clear();
        currentPair = data.at(kk);
        currentList = currentPair.first;
        for (unsigned int k = 0; k < currentList.size(); ++k) {
            for (unsigned int i = 0; i < currentList.at(k).size(); ++i) { //Здесь только те тела, которые +1 или -1.
                zoneVector.push_back(currentList.at(k).at(i).first);
                bodyVector.push_back(currentList.at(k).at(i).second);
            }
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


} //end of namespace geant4toshield

Geant4ShieldData convertGeant4ToShield(G4RunManager *runManager) {
    using namespace geant4toshield;
    G4VUserDetectorConstruction *detector = (G4VUserDetectorConstruction *) runManager->GetUserDetectorConstruction(); //from const to non-const pointer
    G4VPhysicalVolume *physicalWorld = detector->Construct();
    zoneData zoneWorld = convertBodyRecursive(physicalWorld);
    zoneWorld.push_back(addOuterVacuum(physicalWorld));
    return zoneDataToShield(zoneWorld);
} //TODO сделать нормально через map: рекурсивно для всех дочерних получить параметры и внести в map, после этого для получения тела уже брать оттуда.
