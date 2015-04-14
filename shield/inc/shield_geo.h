// File:    shield_geo.h
// Authors: Alexander Timofeev    02/2012
//          Dmitry    Sosnov      02/2015

#ifndef SHIELD_H
#error ERROR: Don`t include shield_geo.h, include shield.h instead.
#endif

#ifndef SHIELD_GEO_H
#define SHIELD_GEO_H

// Geometry callbacks

typedef float (*GEODISTCALLBACK)(float x, float y, float z, 
                                 float vx, float vy, float vz);

typedef int (*GEONEXTCALLBACK)(float x, float y, float z, 
                               float vx, float vy, float vz);

struct SGeoBody {
    int type; //Warning! Type is wrote before parameters in bodydb
    float parameters[36];
};
struct SGeoZone {
    int countELements;
    int mediaNumber;
    int definition[10000]; //Current count of elements at fortran ZONEDB is 5000, but we hawe number of body and type of zone separately.
};

//typedef int (*GEOMEDIUMCALLBACK)(

SHIELD_API void shield_geo_set_dist_callback(GEODISTCALLBACK function);
SHIELD_API void shield_geo_set_next_callback(GEONEXTCALLBACK function);

int shield_add_body(int type, float *parameters);
SHIELD_API int shield_add_zone(int countBodies,  int *zoneParameters, struct SGeoBody *bodies, struct MediumData medium);
SHIELD_API void shield_clean_geometry();
void shield_init_geometry();

#endif
