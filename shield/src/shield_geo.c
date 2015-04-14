#define SHIELD_LIB_INTERNAL
#include "shield.h"
#include "shield_common_blocks.h"
#include <stdio.h>
#include <math.h>
GEODISTCALLBACK geo_dist_callback = 0;
GEONEXTCALLBACK geo_next_callback = 0;

#ifndef gdatap
#define gdatap gdatap_
struct Gdatap_t {
    double x, y, z, cx, cy, cz;
    double xn, yn, zn, cxn, cyn, czn;
    int NZONO, MEDOLD, NBPO, NSO;
    int NZONC, MEDCUR, NBPC, NSCI, NSCO;
    int PINSFL, PARTZD;
} extern gdatap;
#endif

/*--------------------------------------------------------------+*/
/*                     1   2   3   4   5   6   7   8   9  |  10 |*/
/* Number at struct    0   1   2   3   4   5   6   7   8  |  9  |*/
/* Body               SPH WED ARB BOX RPP RCC REC TRC ELL | END |*/
/* N Parameters        6  30  36  36  18  21  22  22  17  |  0  |*/
/* N Surfaces          1   5   6   6   6   3   3   3   1  |  0  |*/
/*--------------------------------------------------------------+*/
/*--------------------------------------------------------------+*/
/* Surface type:      Equation:                 N Parameters:   |*/
/*             1      X-C=0                                 1   |*/
/*             2      Y-C=0                                 1   |*/
/*             3      Z-C=0                                 1   |*/
/*             4      A*x+B*y+C*z-D=0                       4   |*/
/*             5      (x-x0)**2+(y-y0)**2+(z-z0)**2-R**2=0  4   |*/
/*             6      x**2+y**2-R**2=0                      1   |*/
/*             7      x**2/A2+y**2/B2-1=0                   2   |*/
/*             8      x**2+y**2-(z-A1)**2/B2=0              2   |*/
/*             9      x**2/A2+y**2/B2+z**2/C2-1=0           3   |*/
/*--------------------------------------------------------------+*/
struct {
    int count; //count of types of bodies
    int countOfParameters[9];
    char namesOfBodies[9][3];
} ShieldBodies = {9, {6, 30, 36, 36, 18, 21, 22, 22, 17}, {"SPH", "WED", "ARB", "BOX", "RPP", "RCC", "REC", "TRC", "ELL"}};
int ShieldBodiesMaxCountOfParameters() {
    int i, max = 0;
    for (i = 0; i < ShieldBodies.count; i++)
        if (ShieldBodies.countOfParameters[i] > max) {
            max = ShieldBodies.countOfParameters[i];
        }
    return max;
}

float DIGORT = 0.000001;

struct {
    int countBodies;
    struct SGeoBody bodies [1000];
    int countZones;
    struct SGeoZone zones [1000];
} SGeometry = {0, {}, 0, {}};

void shield_geo_set_dist_callback(GEODISTCALLBACK function) {
    geo_dist_callback = function;
}

void shield_geo_set_next_callback(GEONEXTCALLBACK function) {
    geo_next_callback = function;
}

void gcurzl(int *zone_in, float *dist_out) {
    if (geo_dist_callback != 0) {
        *dist_out = geo_dist_callback(gdatap.x, gdatap.y, gdatap.z, gdatap.cx, gdatap.cy, gdatap.cz);
    } else {
        // should never happen
        printf("FATAL ERROR: GEO DIST CALLBACK NOT SET\n");
//        exit(1);
    }
}

void gnextz(int *next_out) {
    // determinition of next zone
    if (geo_next_callback != 0) {
        *next_out = geo_next_callback(gdatap.x, gdatap.y, gdatap.z, gdatap.cx, gdatap.cy, gdatap.cz);
    } else {
        // should never happen
        printf("FATAL ERROR: GEO NEXT CALLBACK NOT SET\n");
//        exit(1);
    }
}

void gzmed(int *zone) {
    // Determinition of zone medium
}

void assigmentToVec3(float *vec, float x1, float x2, float x3) {
    vec[0] = x1;
    vec[1] = x2;
    vec[2] = x3;
}
void assigmentToVec4(float *vec, float x1, float x2, float x3, float x4) {
    vec[0] = x1;
    vec[1] = x2;
    vec[2] = x3;
    vec[3] = x4;
}

//gemca.f:1602
void geqpln(float x, float y, float z, float v1, float v2, float v3, float *out) {
    float d = v1 * x + v2 * y + v3 * z;
    if (d == 0) {
        assigmentToVec4(out, v1, v2, v3, 0);
    } else {
        assigmentToVec4(out, v1 / d, v2 / d, v3 / d, 1);
    }
//     *out = (d == 0) ? {v1,v2,v3,0} : {v1/d,v2/d,v3/d,1};
}

//gemca.f:1579
void gcosax(float x0, float y0, float z0, float x, float y, float z,
            float *out, float *r, float *r2) {
    float *cosx = &out[0];
    float *cosy = &out[1];
    float *cosz = &out[2];
    float *delta = &out[3];
    *r2 = x * x + y * y + z * z;
    *r = sqrt(*r2);
    *cosx = x / (*r);
    *cosy = y / (*r);
    *cosz = z / (*r);
    *delta = x0 * (*cosx) + y0 * (*cosy) + z0 * (*cosz);
}

//gemca.f:1626
//bool nv: if nv.eq.1, then true, if nv.eq.2 then false
//gvec90_1 is gvec90 with nv==1;
//gvec90_2 is gvec90 with nv==2;
void gvec90_1(float x1, float y1, float z1,
              float x2, float y2, float z2, float *out) {
    float a, b, c;
    if (x1 == 0 && x2 == 0) {
        assigmentToVec3(out, 1, 0, 0);
    } else if (y1 == 0 && y2 == 0) {
        assigmentToVec3(out, 0, 1, 0);
    } else if (z1 == 0 && z2 == 0) {
        assigmentToVec3(out, 0, 0, 1);
    } else if (x1 == 0 && y1 == 0) {
        assigmentToVec3(out, -y2, x2, 0);
    } else if (x1 == 0 && z1 == 0) {
        assigmentToVec3(out, -z2, 0, x2);
    } else if (y1 == 0 && z1 == 0) {
        assigmentToVec3(out, 0, -z2, y2);
    } else if (x2 == 0 && y2 == 0) {
        assigmentToVec3(out, -y1, x1, 0);
    } else if (x2 == 0 && z2 == 0) {
        assigmentToVec3(out, -z1, 0, x1);
    } else if (y2 == 0 && z2 == 0) {
        assigmentToVec3(out, 0, -z1, y1);
    } else if (x1 != 0 && x2 != 0) {
        if (y1 / x1 == y2 / x2) {
            assigmentToVec3(out, y1, -x1, 0);
        } else {
            float c = 1;
            float b = (z2 * x1 - z1 * x2) / (y1 * x2 - y2 * x1);
            float a = -(b * y1 + z1) / x1;
            assigmentToVec3(out, a, b, c);
        }
    } else if (y1 != 0 && y2 != 0) {
        if (z1 / y1 == z2 / y2) {
            assigmentToVec3(out, 0, -z1, y1);
        } else {
            float a = 1;
            float c = (x2 * y1 - x1 * y2) / (z1 * y2 - z2 * y1);
            float b = -(c * z1 + x1) / y1;
            assigmentToVec3(out, a, b, c);
        }
    } else if (z1 != 0 && z2 != 0) {
        if (x1 / z1 == x2 / z2) {
            assigmentToVec3(out, -z1, 0, x1);
        } else {
            float b = 1;
            float a = (y2 * z1 - y1 * z2) / (x1 * z2 - x2 * z1);
            float c = -(a * x1 + y1) / z1;
            assigmentToVec3(out, a, b, c);
        }
    }
}
void gvec90_2(float x1, float y1, float z1,
              float *x2, float *y2, float *z2,
              float *out1, float *out2, float *out3) {
    float out[3];
    if (x1 == 0) {
        *x2 = 1;
    } else {
        *x2 = -z1;
    }
    if (y1 == 0) {
        *y2 = 1;
    } else {
        *y2 = 0;
    }
    if (z1 == 0) {
        *z2 = 1;
    } else {
        *z2 = x1;
    }
    gvec90_1(x1, y1, z1, *x2, *y2, *z2, &out[0]);
    *out1 = out[0];
    *out2 = out[1];
    *out3 = out[2];
}
//gemca.f:857
float gssign(int ksurf, float *bodydb, float x, float y, float z) {
    float ssign, trues;
    switch (ksurf) {
        case 1: {
            float c = bodydb[0];
            trues = bodydb[1];
            ssign = x - c;
            break;
        }
        case 2: {
            float c = bodydb[0];
            trues = bodydb[1];
            ssign = y - c;
            break;
        }
        case 3: {
            float c = bodydb[0];
            trues = bodydb[1];
            ssign = z - c;
            break;
        }
        case 4: {
            float a = bodydb[0];
            float b = bodydb[1];
            float c = bodydb[2];
            float d = bodydb[3];
            trues = bodydb[4];
            ssign = a * x + b * y + c * z - d;
            break;
        }
        case 5: {
            float x0 = x - bodydb[0];
            float y0 = y - bodydb[1];
            float z0 = z - bodydb[2];
            float r2 = bodydb[3];
            trues = bodydb[4];
            ssign = (x0 * x0 + y0 * y0 + z0 * z0) / r2 - 1;
            break;
        }
        case 6: {
            float r2 = bodydb[0];
            trues = bodydb[1];
            ssign = x * x + y * y - r2;
            break;
        }
        case 7: {
            float a2 = bodydb[0];
            float b2 = bodydb[1];
            trues = bodydb[2];
            ssign = x * x / a2 + y * y / b2 - 1;
            break;
        }
        case 8: {
            float a1 = bodydb[0];
            float b2 = bodydb[1];
            trues = bodydb[2];
            ssign = x * x + y * y - (z - a1) * (z - a1) / b2;
            break;
        }
        case 9: {
            float a2 = bodydb[0];
            float b2 = bodydb[1];
            float c2 = bodydb[2];
            trues = bodydb[3];
            ssign = x * x / a2 + y * y / b2 + z * z / c2 - 1;
            break;
        }
        case 10:
            printf("EXIT from GSSIGN : no such surface !");
            break;
    }
    //TODO Проверить выходные значения
    return (ssign * trues);
}

int testOrtoganality(int body, float *point1, float *point2) {
    int i;
    float dort = point1[0] * point2[0] + point1[1] * point2[1] + point1[2] * point2[2];
    if (fabs(dort) >= DIGORT) {
        printf("%.3s: error: not orthogonal, abs=%f\n", ShieldBodies.namesOfBodies[body], fabs(dort));
        return 1;
    }
    return 0;
}
//Parameters is coordinates(x,y,z) and radius
struct SGeoBody createSPH(float *parameters) {
    float x = parameters[0];
    float y = parameters[1];
    float z = parameters[2];
    float r = parameters[3];
    struct SGeoBody body = {0, {5, x, y, z, r * 2, -1}};
    return body;
}
//Parameters is coords of 4 points
struct SGeoBody createWED(float *parameters) {
    struct SGeoBody body = {1, {}};
    float tmpxyz [3];
    int i;
    for (i = 0; i < 2; i++) {
        if (testOrtoganality(1, &parameters[3], &parameters[3 * (i + 2)])) {
            return;
        }
    }
    float x = parameters[0] + 0.25 * (parameters[3] + parameters[6] + parameters[9]);
    float y = parameters[1] + 0.25 * (parameters[4] + parameters[7] + parameters[10]);
    float z = parameters[2] + 0.25 * (parameters[5] + parameters[8] + parameters[11]);

    geqpln(parameters[0], parameters[1], parameters[2],
           parameters[3], parameters[4], parameters[5],
           &body.parameters[1]);
    geqpln(parameters[0] + parameters[3], parameters[1] + parameters[4], parameters[2] + parameters[5],
           parameters[3], parameters[4], parameters[5],
           &body.parameters[7]);
    gvec90_1(parameters[3], parameters[4], parameters[5], parameters[6], parameters[7], parameters[8], &tmpxyz[0]);
    geqpln(parameters[0], parameters[1], parameters[2],
           tmpxyz[0], tmpxyz[1], tmpxyz[2],
           &body.parameters[13]);
    gvec90_1(parameters[3], parameters[4], parameters[5], parameters[9], parameters[10], parameters[11], &tmpxyz[0]);
    geqpln(parameters[0], parameters[1], parameters[2],
           tmpxyz[0], tmpxyz[1], tmpxyz[2],
           &body.parameters[19]);
    gvec90_1(parameters[3], parameters[4], parameters[5],
             parameters[9] - parameters[6], parameters[10] - parameters[7], parameters[11] - parameters[8],
             &tmpxyz[0]);
    geqpln(parameters[0] + parameters[6], parameters[1] + parameters[7], parameters[2] + parameters[8],
           tmpxyz[0], tmpxyz[1], tmpxyz[2],
           &body.parameters[25]);

    for (i = 0; i < 5; i++) {
        body.parameters[i * 6] = 4;
        body.parameters[i * 6 + 5] = 1;
        body.parameters[i * 6 + 5] = gssign(4, &body.parameters[i * 6 + 1], x, y, z);
    }
    return body;
}

//Positions - list of numbers for (1-8)
struct SGeoBody createARBWithPositions(float *parameters, int *positions) {
    struct SGeoBody body = {2, {}};
    int i;
    float x = 0, y = 0, z = 0;
    float tmpxyz [3];
    for (i = 0; i < 8; i++) {
        x += parameters[3 * i];
        y += parameters[3 * i + 1];
        z += parameters[3 * i + 2];
    }
    x /= 8; y /= 8; z /= 8;
    for (i = 0; i < 6; i++) {
        body.parameters[i * 6] = 4;
        int ip1 = (positions[i * 4] - 1) * 3;
        int ip2 = (positions[i * 4 + 1] - 1) * 3;
        int ip3 = (positions[i * 4 + 2] - 1) * 3;
        int ip4 = (positions[i * 4 + 3] - 1) * 3;
        gvec90_1(parameters[ip1] - parameters[ip2], parameters[ip1 + 1] - parameters[ip2 + 1], parameters[ip1 + 2] - parameters[ip2 + 2],
                 parameters[ip1] - parameters[ip3], parameters[ip1 + 1] - parameters[ip3 + 1], parameters[ip1 + 2] - parameters[ip3 + 2],
                 &tmpxyz[0]);
        geqpln(parameters[ip1], parameters[ip1 + 1], parameters[ip1 + 1],
               tmpxyz[0], tmpxyz[1], tmpxyz[2],
               &body.parameters[i * 6 + 1]);
        body.parameters[i * 6 + 5] = 1;
        body.parameters[i * 6 + 5] = gssign(4, &body.parameters[i * 6 + 1], x, y, z);
        float checks = gssign(4, &body.parameters[i * 6 + 1], parameters[ip4], parameters[ip4 + 1], parameters[ip4 + 2]);
        if (abs(checks) >= DIGORT) { //In gemca was DIGMIN
            printf("EXIT from GEOINI: Invalid coordinates !");
            return;
        }
    }
    return body;
}
struct SGeoBody createARB(float *parameters) {
    int positions[24] = {
        1, 2, 3, 4,
        5, 1, 4, 8,
        8, 4, 3, 7,
        7, 3, 2, 6,
        6, 2, 1, 5,
        6, 5, 8, 7
    };
    return createARBWithPositions(parameters, &positions[0]);
}
struct SGeoBody createBOX(float *parameters) {
    struct SGeoBody body = {3, {}};
    float tmpxyz [3];
    int i;
    for (i = 0; i < 2; i++) {
        if (testOrtoganality(3, &parameters[3], &parameters[3 * (i + 2)])) {
            return;
        }
    }
    if (testOrtoganality(3, &parameters[6], &parameters[9])) {
        return;
    }
    float x = parameters[0] + 0.5 * (parameters[3] + parameters[6] + parameters[9]);
    float y = parameters[1] + 0.5 * (parameters[4] + parameters[7] + parameters[10]);
    float z = parameters[2] + 0.5 * (parameters[5] + parameters[8] + parameters[11]);

    int k;
    for (i = 0; i < 3; i++) {
        k = i * 3 + 3;
        geqpln(parameters[0], parameters[1], parameters[2],
               parameters[k], parameters[k + 1], parameters[k + 2],
               &body.parameters[i * 12 + 1]);
        geqpln(parameters[0] + parameters[k], parameters[1] + parameters[k + 1], parameters[2] + parameters[k + 2],
               parameters[k], parameters[k + 1], parameters[k + 2],
               &body.parameters[i * 12 + 7]);
    }
    for (i = 0; i < 6; i++) {
        body.parameters[i * 6] = 4;
        body.parameters[i * 6 + 5] = 1;
        body.parameters[i * 6 + 5] = gssign(4, &body.parameters[i * 6 + 1], x, y, z);
    }
    return body;
}
struct SGeoBody createRPP(float *parameters) {
    struct SGeoBody body = {4, {}};
    if (parameters[0] >= parameters[1] ||
            parameters[2] >= parameters[3] ||
            parameters[4] >= parameters[5]) {
        printf("%s: X1 must be less than X2", ShieldBodies.namesOfBodies[4]);
        printf("%s: Y1 must be less than Y2", ShieldBodies.namesOfBodies[4]);
        printf("%s: Z1 must be less than Z2", ShieldBodies.namesOfBodies[4]);
        return;
    }
    int i;
    for (i = 0; i < 3; i++) {
        body.parameters[i * 6] = i + 1;
        body.parameters[i * 6 + 1] = parameters[i * 2];
        body.parameters[i * 6 + 2] = 1;
        body.parameters[i * 6 + 3] = i + 1;
        body.parameters[i * 6 + 4] = parameters[i * 2 + 1];
        body.parameters[i * 6 + 5] = -1;
    }
    return body;
}
struct SGeoBody createRCC(float *parameters) {
    struct SGeoBody body = {5, {}};
    float a1, b1, c1, a2, b2, c2, r, r2;
    gvec90_2(parameters[3], parameters[4], parameters[5],
             &a1, &b1, &c1, &a2, &b2, &c2);
    gcosax(parameters[0], parameters[1], parameters[2],
           a1, b1, c1,
           &body.parameters[0], &r, &r2);
    gcosax(parameters[0], parameters[1], parameters[2],
           a2, b2, c2,
           &body.parameters[4], &r, &r2);
    gcosax(parameters[0], parameters[1], parameters[2],
           parameters[3], parameters[4], parameters[5],
           &body.parameters[8], &body.parameters[16], &r2);
//             &body.parameters[8],r,r2);
    body.parameters[12] = 3;
    body.parameters[13] = 0;
    body.parameters[14] = 1;
    body.parameters[15] = 3;
//     body.parameters[16] = r;
    body.parameters[17] = -1;
    body.parameters[18] = 6;
    body.parameters[19] = parameters[6] * parameters[6];
    body.parameters[20] = -1;
    return body;
}
struct SGeoBody createREC(float *parameters) {
    struct SGeoBody body = {6, {}};
    float r11, r12, r21, r22, h1, h2;
    int i;
    for (i = 0; i < 2; i++) {
        if (testOrtoganality(6, &parameters[3], &parameters[3 * (i + 2)])) {
            return;
        }
    }
    gcosax(parameters[0], parameters[1], parameters[2],
           parameters[6], parameters[7], parameters[8],
           &body.parameters[0], &body.parameters[19], &r12);
    gcosax(parameters[0], parameters[1], parameters[2],
           parameters[9], parameters[10], parameters[11],
           &body.parameters[4], &body.parameters[20], &r22);
    gcosax(parameters[0], parameters[1], parameters[2],
           parameters[3], parameters[4], parameters[5],
           &body.parameters[8], &body.parameters[16], &h2);
    body.parameters[12] = 3;
    body.parameters[13] = 0;
    body.parameters[14] = 1;
    body.parameters[15] = 3;
//     body.parameters[16] = h1;
    body.parameters[17] = -1;
    body.parameters[18] = 7;
//     body.parameters[19] = r12;
//     body.parameters[20] = r22;
    body.parameters[21] = -1;
    return body;
}
struct SGeoBody createTRC(float *parameters) {
    struct SGeoBody body = {7, {}};
    float r, r2, h, h2;
    float a1, b1, c1, a2, b2, c2;
    if (parameters[7] >= parameters[6]) {
        printf("%s: R2 must be less than R1", ShieldBodies.namesOfBodies[7]);
        return;
    }
    gvec90_2(parameters[3], parameters[4], parameters[5],
             &a1, &b1, &c1, &a2, &b2, &c2);
    gcosax(parameters[0], parameters[1], parameters[2],
           a1, b1, c1,
           &body.parameters[0], &r, &r2);
    gcosax(parameters[0], parameters[1], parameters[2],
           a2, b2, c2,
           &body.parameters[4], &r, &r2);
    gcosax(parameters[0], parameters[1], parameters[2],
           parameters[3], parameters[4], parameters[5],
           &body.parameters[8], &h, &h2);
    body.parameters[12] = 3;
    body.parameters[13] = 0;
    body.parameters[14] = 1;
    body.parameters[15] = 3;
    body.parameters[16] = h;
    body.parameters[17] = -1;
    body.parameters[18] = 8;
    body.parameters[19] = h / (1 - parameters[7] - parameters[6]);
    body.parameters[20] = (body.parameters[19] * body.parameters[19]) / (parameters[6] * parameters[6]);
    body.parameters[21] = -1;
    return body;
}
struct SGeoBody createELL(float *parameters) {
    struct SGeoBody body = {7, {}};
    float r;
    int i;
    for (i = 0; i < 2; i++) {
        if (testOrtoganality(7, &parameters[3], &parameters[3 * (i + 2)])) {
            return;
        }
    }
    if (testOrtoganality(7, &parameters[6], &parameters[9])) {
        return;
    }
    for (i = 0; i < 3; i++) {
        gcosax(parameters[0], parameters[1], parameters[2],
               parameters[3 * (i + 1)], parameters[3 * (i + 1) + 1], parameters[3 * (i + 1) + 2],
               &body.parameters[i * 4], &r, &body.parameters[13 + i]);
    }
    body.parameters[12] = 9;
    body.parameters[16] = -1;
    return body;
}
struct SGeoBody createBody(int type, float *parameters) {
    switch (type) {
        case 0:
            return createSPH(parameters);
            break;
        case 1:
            return createWED(parameters);
            break;
        case 2:
            return createARB(parameters);
            break;
        case 3:
            return createBOX(parameters);
            break;
        case 4:
            return createRPP(parameters);
            break;
        case 5:
            return createRCC(parameters);
            break;
        case 6:
            return createREC(parameters);
            break;
        case 7:
            return createTRC(parameters);
            break;
        case 8:
            return createELL(parameters);
            break;
        default:
            printf("Error in number of body");
            struct SGeoBody defaultBody = { -1, {}};
            return defaultBody;
            break;
    }
}

//1 is True and 0 is false
int isBodiesEqual(struct SGeoBody body1, struct SGeoBody body2) {
    int i;
    if (body1.type != body2.type) {
        return 0;
    }
    for (i = 0; i < ShieldBodies.countOfParameters[body1.type]; i++)
        if (body1.parameters[i] != body2.parameters[i]) {
            return 0;
        }
    return 1;
}

int findBody(struct SGeoBody body) {
    int i;
    for (i = 0; i < SGeometry.countBodies; i++) {
        if (isBodiesEqual(body, SGeometry.bodies[i])) {
            return i;
        }
    }
    return -1;
}
int shield_add_body(int type, float *parameters) {
    struct SGeoBody body = createBody(type, parameters);
    if (body.type == -1) {
        return -1;
    }
    int bodynum = findBody(body);
    if (bodynum != -1) {
        return bodynum;
    }
    if(SGeometry.countBodies>=1000){ //Parameter MAXB at fortran code
        printf("shield_init_geometry ERROR: Count of bodies large then array for bodies.\n");
        return;
    }
    SGeometry.bodies[SGeometry.countBodies] = body;
    return SGeometry.countBodies++; //returned SGeometry.countBodies before increment
}

//if testMedia == 0, then only geometry testing
int isZonesEqual(int testMedia, struct SGeoZone zone1, struct SGeoZone zone2) {
    int i;
    if (zone1.countELements != zone2.countELements) {
        return 0;
    }
    if (testMedia && zone1.mediaNumber != zone2.mediaNumber) {
        return 0;
    }
    for (i = 0; i < zone1.countELements * 2; i++)
        if (zone1.definition[i] != zone2.definition[i]) {
            return 0;
        }
    return 1;
}

//Return value is boolean.
int shield_add_zone(int countBodies,  int *zoneParameters, struct SGeoBody *bodies, struct MediumData medium) {
    int i, c;
    int bodynum;
    struct SGeoZone zone = {0, 0, {}};
    for (i = 0, c = 0 ; c < countBodies; i++) {
        zone.definition[zone.countELements * 2] = zoneParameters[i];
        switch (zoneParameters[i]) {
            case 0:
                break;
            case 1:
            case -1:
                bodynum = shield_add_body(bodies[c].type, bodies[c].parameters);
                if (bodynum == -1) {
                    printf("shield_init_geometry ERROR: Cannot add body\n");
                    return 1;
                }
                zone.definition[zone.countELements * 2 + 1] = bodynum;
                c++;
                break;
            default:
                printf("shield_init_geometry ERROR: Uncorrect array for zone\n");
                return 1;
                break;
        }
        zone.countELements++;
        if(zone.countELements>5000){ // Parameter, count of elements at array at struct
            printf("shield_init_geometry ERROR: Count of elements at zone defenition more then capacity\n");
            return 1;
        }
    }
    int zonenum = -1;
    for (i = 0; i < SGeometry.countZones; i++) {
        if (isZonesEqual(0, zone, SGeometry.zones[i])) {
            zonenum = i;
            break;
        }
    }
    if (zonenum != -1) {
        printf("shield_init_geometry WARNING: Doublicate of zone, skipping\n");
        return 0;
    }
    if (medium.nType == 0 || medium.nType == 1000) {
        zone.mediaNumber = medium.nType-1;
    } else {
        zone.mediaNumber = shield_add_medium(medium);
    }

    if(SGeometry.countZones>=1000){ //Parameter MAXZ at fortran code
        printf("shield_init_geometry ERROR: Count of zones large then array for zones.\n");
        return 1;
    }
    SGeometry.zones[SGeometry.countZones] = zone;
    SGeometry.countZones++;
//     return SGeometry.countZones++; //returned SGeometry.countBodies before increment
    return 0;
}

void shield_init_geometry() {
    int i, j;
//     printf("\tCurrent count of elements: %i, %i\n", SGeometry.countBodies, SGeometry.countZones);
    gdata0.NUMBOD = SGeometry.countBodies;
    gdata0.NUMZON = SGeometry.countZones;
    int curPos = 0;
    struct SGeoBody body;
    struct SGeoZone zone;
    for (i = 0; i < SGeometry.countBodies; i++) {
        body = SGeometry.bodies[i];
        gdata2.BODYDA[i] = curPos + 1;
        gdata2.BODYDB[curPos] = body.type + 1;
        curPos++;
        for (j = 0; j < ShieldBodies.countOfParameters[body.type]; j++) {
            gdata2.BODYDB[curPos + j] = body.parameters[j];
        }
        curPos += ShieldBodies.countOfParameters[body.type];
    }
    curPos = 0;
    for (i = 0; i < SGeometry.countZones; i++) {
        zone = SGeometry.zones[i];
        if(curPos + zone.countELements + 3 > 5000){
            printf("shield_init_geometry ERROR: Too large zones.\n");
            return;
        }
        gdata3.ZONEDA[i] = curPos + 1;
        gdata3.ZONEDB[curPos] = zone.countELements;
        gdata3.ZONEDB[curPos + 1] = zone.mediaNumber + 1;
        int typeOfElements = 1;
        for (j = 0; j < zone.countELements; j++) {
            if (zone.definition[j * 2] == 0) {
                typeOfElements = 0;
            } else if (typeOfElements != 0 && zone.definition[j * 2] == -1) { //in gemca:1515
                typeOfElements = -1;
            }
            gdata3.ZONEDB[curPos + 3 + j] = zone.definition[j * 2] * (zone.definition[j * 2 + 1] + 1);
        }
        gdata3.ZONEDB[curPos + 2] = typeOfElements;
        curPos += zone.countELements + 3;
    }
}
void shield_clean_geometry(){
    SGeometry.countBodies=0;
    SGeometry.countZones=0;
    shield_clean_medium();
}
