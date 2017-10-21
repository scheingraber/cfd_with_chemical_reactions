#ifndef __REACTION_H_
#define __REACTION_H_

#include <math.h>
#include <algorithm>
#include "Parameters.h"
#include <float.h>
#include <assert.h>

// This function has to produce a time step for every component, at every
// point, given all the reactions, that cannot produce a future negative value
// for the concentration.
double reaction_max_dt(
        double ***C,
        double **T,
        int **Flag,
        const Parameters & params,
        std::vector<double> & rates
        );

void compute_reaction(
        double ***C,
        double **T,
        int **Flag,
        double dt,
        const Parameters & params,
        std::vector<double> & rates
        );

#endif
