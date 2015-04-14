// File:    shield_tree.h
// Authors: Alexander Timofeev 10/02/2012
//          Dmitry    Sosnov      02/2015

#ifndef SHIELD_TREE_H
#define SHIELD_TREE_H

#ifndef SHIELD_H
#error ERROR: Don`t include shield_tree.h header, include shield.h instead
#endif

// base structures of TREE data
// should be used to derive classes to import data to ROOT

#define EVENT_PRIMARY           0
#define EVENT_ELASTIC           1
#define EVENT_INELASTIC         2
#define EVENT_FLYOUT            3
#define EVENT_ABSORPTION        4
#define EVENT_DECAY             5
#define EVENT_DECAYREST         6

struct shield_tree_node {
    int event;              // event type (see above 0..6 definitions)
    int jpart;              // particle type (see hadgen documentation)
    int pdg;                // PDG code of particle or nucleus
    float weight;           // particle weight [GeV/c^2]
    int nbranch;            // branch number
    float xi, yi, zi, ti;   // initial position of branch [m]
    float xf, yf, zf, tf;   // point of interest in event (star or elastic scattering)
    float cost, sinf, cosf; // direction        
    float xz, yz, zz, tz;   // final point -- equal to xf.. in all cases except elastic star
    float px_z, py_z, pz_z; // momentum (computed from direction) [GeV/c]
    int iapr, izpr;         // A and Z of nucleus
    int izsc, iasc;
    int istate, nbrnew;     // final state, next branch number 
                            //              if applicable

};

// callback function to external storage
// NOTICE: this callback should NOT store the pointer
// the pointer will be freed just after callback call
typedef void (*TREECALLBACK)(struct shield_tree_node *);

SHIELD_API void shield_set_tree_callback(TREECALLBACK func);

#endif
