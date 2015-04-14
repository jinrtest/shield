// File:    shield_tree.h
// Authors: Alexander Timofeev 09.01.2012
// based on original FORTRAN code

#define SHIELD_LIB_INTERNAL
#include "shield.h"
#include "hadgen.h"
#include "shield_common_blocks.h"
#include <math.h>
#include <stdio.h>

// HADGEN internal functions

extern float _hadgen_get_particle_weight(int type);

// Internal data
// Tree callback function

TREECALLBACK tree_callback = 0; 
// WARNING: should be changed in future for parallelism
// (should be included into running context)

void shield_tree_elastic(struct shield_tree_node *node,
        int j_start, int j_end, 
        int n_p_cas, int n_p_tot);
void shield_tree_inelastic(struct shield_tree_node *node,
        int j_start, int j_end, 
        int n_p_cas, int n_p_tot);
void shield_tree_decay(struct shield_tree_node *node,
        int j_start, int j_end, 
        int n_p_cas, int n_p_tot);
void shield_tree_decayrest(struct shield_tree_node *node,
        int j_start, int j_end, 
        int n_p_cas, int n_p_tot);
void shield_tree_primary(struct shield_tree_node *node);
void shield_tree_flyout(struct shield_tree_node *node,
        int j_start, int j_end);
void shield_tree_absorption(struct shield_tree_node *node,
        int j_start, int j_end);

void shield_fill_node(struct shield_tree_node *node)
{
//	printf("jpart %i\n", node->jpart);
    if (node->jpart != 25) {
        node->weight = _hadgen_get_particle_weight(node->jpart);
        node->pdg = hadgen_get_pdg_code(node->jpart);
    } else {
        node->weight = 0.940 * node->iapr; // in [GeV/c^2]
        node->pdg = hadgen_get_pdg_code_nuclei(node->iapr, node->izpr);
    }
//	printf("shield, pdg: %i\n", node->pdg);
    node->ti *= 0.001;
    node->tf *= 0.001;
    node->tz *= 0.001;
    float E = node->tz;
    float P = sqrt(E * (E + 2.0 * node->weight)); // in [GeV/c]
    float sint = sqrt(1.0 - node->cost * node->cost);
    node->px_z = P * sint * node->cosf;
    node->py_z = P * sint * node->sinf;
    node->pz_z = P * node->cost;
}

void shield_set_tree_callback(TREECALLBACK function)
{
    tree_callback = function;
}

void shield_tree_new_branch(int *kstate, // treedec_
    int *j_start, int *j_end, 
    int *n_part_cas, int *n_part_tot)
{
//    printf("treedec_\n");
    if (tree_callback != 0) {
        struct shield_tree_node node;
        // do everything only if callback is set
        switch(*kstate) {
            case 1:         // elastic scattering
                shield_tree_elastic(&node, *j_start-1, *j_end-1, *n_part_cas, *n_part_tot);
                break;
            case 2:         // inelastic star
                shield_tree_inelastic(&node, *j_start-1, *j_end-1, *n_part_cas, *n_part_tot);
                break;
            case 5:         // decay on fly
                shield_tree_decay(&node, *j_start-1, *j_end-1, *n_part_cas, *n_part_tot);
                break;
            case 6:         // decay at rest
                shield_tree_decayrest(&node, *j_start-1, *j_end-1, *n_part_cas, *n_part_tot);
                break;
            case 3:         // fly out
                shield_tree_flyout(&node, *j_start-1, *j_end-1);
                break;
            case 4:         // absorption
                shield_tree_absorption(&node, *j_start-1, *j_end-1);
                break;
            case 0:         // primary particle
                shield_tree_primary(&node);
                break;
            default:
                break;            
        }
        // fill the node before callback
        shield_fill_node(&node);
        tree_callback(&node); // reporting outside...
    }
}

void shield_tree_reset(void)    // treeres_
{
//    printf("treeres_\n");
    // dummy function at the moment
}

void shield_tree_prewrite_data(struct shield_tree_node *node, int j_start)
{
    node->jpart = tree.ITREE[j_start][0]; // current particle type
    node->nbranch = tree.ITREE[j_start][2]; // number of current branch

    // if particle is a nucleus...
    if (node->jpart == 25) {
        // write A and Z of the nucleus
        node->iapr = (int)(tree.TREE[j_start][13]);
        node->izpr = (int)(tree.TREE[j_start][14]);
    } else {
        node->iapr = 0;
        node->izpr = 0;
    }
    // initial position
    node->xi = tree.TREE[j_start][0];
    node->yi = tree.TREE[j_start][1];
    node->zi = tree.TREE[j_start][2];
    // initial energy
    node->ti = tree.TREE[j_start][6];
}

void shield_tree_write_pos(struct shield_tree_node *node, int j_start)
{
    // scattering point
    node->xf = tree.TREE[j_start+1][0];
    node->yf = tree.TREE[j_start+1][1];
    node->zf = tree.TREE[j_start+1][2];
    // energy just before scattering
    node->tf = tree.TREE[j_start+1][6];
    // write final point xz, yz, zz, tz the same as xf, yf, zf, tf
    // this is by default, except elastic case!
    node->xz = node->xf;
    node->yz = node->yf;
    node->zz = node->zf;
    node->tz = node->tf;
    // impulse
    node->cost = tree.TREE[j_start+1][3];
    node->sinf = tree.TREE[j_start+1][4];
    node->cosf = tree.TREE[j_start+1][5];
}

void shield_tree_primary(struct shield_tree_node *node) 
{
    node->event = EVENT_PRIMARY;
    shield_tree_prewrite_data(node, 1);
    node->cost = tree.TREE[1][3];
    node->sinf = tree.TREE[1][4];
    node->cosf = tree.TREE[1][5];

    int k0l, k0s;   // WHAT TODO WITH THIS???
    if ((node->jpart == 10) || (node->jpart == 11)) {
        k0l = tree.TREE[1][9]; // IFIX TODO
        k0s = tree.TREE[1][10];
    }
    
    int istate = node->istate = tree.ITREE[2][6];
    if ((istate > 2) && (istate != 5) && (istate != 6)) {
        shield_tree_write_pos(node, 1);
        // absorption or fly-out
	    switch(istate)
		{
		case 3: node->event = EVENT_FLYOUT;
		    break;
		case 4: node->event = EVENT_ABSORPTION;
		    break;
		default:
		    break;
		}
        // at the point xf, yf, zf, tf
    }
}

void shield_tree_elastic(struct shield_tree_node *node,
        int j_start, int j_end, 
        int n_p_cas, int n_p_tot)
{
    node->event = EVENT_ELASTIC;
    shield_tree_prewrite_data(node, j_start);
    shield_tree_write_pos(node, j_start);

    // scattering on nucleus
    node->izsc = tree.TREE[j_start+2][4];
    node->iasc = tree.TREE[j_start+2][3];

    // final point
    node->xz = trevs.TREVS[0][0];
    node->yz = trevs.TREVS[0][1];
    node->zz = trevs.TREVS[0][2];
    // final energy
    node->tz = trevs.TREVS[0][6];

    node->istate = trevs.ITREVS[0][6]; // final status (next branch type)
    node->nbrnew = trevs.ITREVS[0][2]; // new branch number
}

void shield_tree_import_star(struct shield_tree_node *node, 
        int j_start, int j_end,
        int n_p_cas, int n_p_tot)
{
    // TODO all particles!
    // mat be not necessary
}

void shield_tree_inelastic(struct shield_tree_node *node,
        int j_start, int j_end, 
        int n_p_cas, int n_p_tot)
{
    node->event = EVENT_INELASTIC;
    shield_tree_prewrite_data(node, j_start);
    shield_tree_write_pos(node, j_start);

    // collision with nucleus
    node->izsc = tree.TREE[j_start+2][4];
    node->iasc = tree.TREE[j_start+2][3];

    shield_tree_import_star(node, j_start, j_end, n_p_cas, n_p_tot);    
}

void shield_tree_decay(struct shield_tree_node *node,
        int j_start, int j_end, 
        int n_p_cas, int n_p_tot)
{
    node->event = EVENT_DECAY;
    shield_tree_prewrite_data(node, j_start);
    shield_tree_write_pos(node, j_start);

    shield_tree_import_star(node, j_start, j_end, n_p_cas, n_p_tot);    
}


void shield_tree_decayrest(struct shield_tree_node *node,
        int j_start, int j_end, 
        int n_p_cas, int n_p_tot)
{
    // do the same thing as for ord. decay
    shield_tree_decay(node, j_start, j_end, n_p_cas, n_p_tot);
    node->event = EVENT_DECAYREST;
}

void shield_tree_flyout(struct shield_tree_node *node,
        int j_start, int j_end)
{
    // NOT TESTED, TODO
    node->event = EVENT_FLYOUT;
    shield_tree_prewrite_data(node, j_start);
    shield_tree_write_pos(node, j_start);
}

void shield_tree_absorption(struct shield_tree_node *node,
        int j_start, int j_end)
{
    // NOT TESTED, TODO
    node->event = EVENT_ABSORPTION;
    shield_tree_prewrite_data(node, j_start);
    shield_tree_write_pos(node, j_start);
}
