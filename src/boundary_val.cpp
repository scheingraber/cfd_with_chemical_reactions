#include "boundary_val.h"
#include "boundary_conditions.h"
#include "Parameters.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

/**
 * The boundary values of the problem are set.
 */
void domain_boundary_values(
  const Parameters &parameters,
  double **U,
  double **V,
  double **T,
  double ***C
) {
    // Seems to me like the less expensive way to do is just to write a lot
    // just to avoid a lot of comparisons, so let's begin

    // All these loops have both a component for velocity and temperature


    // Left wall

    if (parameters.wlvp == boundary_condition["no-slip"]){
        for (int j = 1; j <= parameters.jmax; j++){
            U[0][j] = 0;
            V[0][j] = -V[1][j];
        }
    }

    else if (parameters.wlvp == boundary_condition["free-slip"]) {
        for (int j = 1; j <= parameters.jmax; j++){
            U[0][j] = 0;
            V[0][j] = V[1][j];
        }
    }

    else if (parameters.wlvp == boundary_condition["outflow"]) {
        for (int j = 1; j <= parameters.jmax; j++){
            U[0][j] = U[1][j];
            V[0][j] = V[1][j];
        }
    }


    if ( parameters.wlt == boundary_condition["dirichlet"] ){ // Fixed temperature
        for (int j = 1; j <= parameters.jmax; j++){
            T[0][j] = 2 * parameters.tl - T[1][j];
        }
    }
    else if ( parameters.wlt == boundary_condition["neumann"] ){ // Uniform Neumann condition
        for (int j = 1; j <= parameters.jmax; j++){
            T[0][j] = T[1][j] - parameters.dx * parameters.tl;
        }
    }


    // Right wall

    if (parameters.wrvp == boundary_condition["no-slip"]){
        for (int j = 1; j <= parameters.jmax; j++){
            U[parameters.imax][j] = 0;
            V[parameters.imax+1][j] = -V[parameters.imax][j];
        }
    }

    else if(parameters.wrvp == boundary_condition["free-slip"]){
        for (int j = 1; j <= parameters.jmax; j++){
            U[parameters.imax][j] = 0;
            V[parameters.imax+1][j] = V[parameters.imax][j];
        }
    }

    else if (parameters.wrvp == boundary_condition["outflow"]){
        for (int j = 1; j <= parameters.jmax; j++){
            U[parameters.imax][j] = U[parameters.imax-1][j];
            V[parameters.imax+1][j] = V[parameters.imax][j];
        }
    }

    if ( parameters.wrt == boundary_condition["dirichlet"] ){ // Fixed temperature
        for (int j = 1; j <= parameters.jmax; j++){
            T[parameters.imax + 1][j] = 2 * parameters.tr - T[parameters.imax][j];
        }
    }
    else if ( parameters.wrt == boundary_condition["neumann"] ){ // Uniform Neumann condition
        for (int j = 1; j <= parameters.jmax; j++){
            T[parameters.imax + 1][j] = T[parameters.imax][j] - parameters.dx * parameters.tr;
        }
    }


    // Upper wall

    if (parameters.wtvp == boundary_condition["no-slip"]) {
        for (int i = 1; i <= parameters.imax; i++){
            U[i][parameters.jmax+1] = -U[i][parameters.jmax];
            V[i][parameters.jmax] = 0;
        }
    }

    else if(parameters.wtvp == boundary_condition["free-slip"]) {
        for (int i = 1; i <= parameters.imax; i++){
            U[i][parameters.jmax+1] = U[i][parameters.jmax];
            V[i][parameters.jmax] = 0;
        }
    }

    else if (parameters.wtvp == boundary_condition["outflow"]) {
        for (int i = 1; i <= parameters.imax; i++){
            U[i][parameters.jmax+1] = U[i][parameters.jmax];
            V[i][parameters.jmax] = V[i][parameters.jmax-1];
        }
    }

    if ( parameters.wtt == boundary_condition["dirichlet"] ){ // Fixed temperature
        for (int i = 1; i <= parameters.imax; i++){
            T[i][parameters.jmax + 1] = 2 * parameters.tt - T[i][parameters.jmax];
        }
    }
    else if ( parameters.wtt == boundary_condition["neumann"] ){ // Uniform Neumann condition
        for (int i = 1; i <= parameters.imax; i++){
            T[i][parameters.jmax + 1] = T[i][parameters.jmax] - parameters.dy * parameters.tt;
        }
    }


    // Bottom wall

    if (parameters.wbvp == boundary_condition["no-slip"]) {
        for (int i = 1; i <= parameters.imax; i++){
            U[i][0] = -U[i][1];
            V[i][0] = 0;
        }
    }

    else if(parameters.wbvp == boundary_condition["free-slip"]) {
        for (int i = 1; i <= parameters.imax; i++){
            U[i][0] = U[i][1];
            V[i][0] = 0;
        }
    }

    else if (parameters.wbvp == boundary_condition["outflow"]) {
        for (int i = 1; i <= parameters.imax; i++){
            U[i][0] = U[i][1];
            V[i][0] = V[i][1];
        }
    }

    if ( parameters.wbt == boundary_condition["dirichlet"] ){ // Fixed temperature
        for (int i = 1; i <= parameters.imax; i++){
            T[i][0] = 2 * parameters.tb - T[i][1];
        }
    }
    else if ( parameters.wbt == boundary_condition["neumann"] ){ // Uniform Neumann condition
        for (int i = 1; i <= parameters.imax; i++){
            T[i][0] = T[i][1] - parameters.dy * parameters.tb;
        }
    }

    // Substances. Taken as zero-flux.
    // Inputs specified in special boundary values according to the problem
    for (unsigned int s = 0 ; s < parameters.nof_substances(); s++){
        for (int i = 1 ; i <= parameters.imax ; i++){
            C[s][i][parameters.jmax + 1] = C[s][i][parameters.jmax];
            C[s][i][0] = C[s][i][1];
        }
        for (int j = 1 ; j <= parameters.jmax ; j++){
            C[s][parameters.imax + 1][j] = C[s][parameters.imax][j];
            C[s][0][j] = C[s][1][j];
        }
    }

}



void spec_boundary_val(
        const char *problem,
        const Parameters & parameters,
        double **U,
        double **V,
        double ***C
) {
    if (!strcmp(problem, "Wire")) { // Flow around a wire
        // loop over the inflow boundary
        for (int j = 1; j <= parameters.jmax; ++j) {

            // Normalization factor for the profile
            double norm_factor;
            if ( parameters.jmax % 2 ){ // Odd number of cells
                norm_factor = ((parameters.jmax / 2) + 0.5) * ((parameters.jmax / 2) + 0.5);
            }
            else{  // Even number of cells
                norm_factor = (parameters.jmax / 2) * (parameters.jmax / 2);
            }

            // Assuming that the velocity is constant and equal to 1
            U[0][j] = 3/2 * (j - 0.5) * ( parameters.jmax - j + 0.5 ) / norm_factor;
            V[0][j] = 0.0;
        }
    }

    if (!strcmp(problem,"Mixing")){
        // Fix parabollic input velocity
        for (int j = 1; j <= parameters.jmax; ++j) {
            double norm_factor;
            if ( parameters.jmax % 2 ){
                norm_factor = ((parameters.jmax / 2) + 0.5) * ((parameters.jmax / 2) + 0.5);
            }
            else{
                norm_factor = (parameters.jmax / 2) * (parameters.jmax / 2);
            }
            U[0][j] = 3/2 * (j - 0.5) * ( parameters.jmax - j + 0.5 ) / norm_factor;
            V[0][j] = 0.0;
        }
        // Then fix constant concentrations at the input
        for (int j = 1; j <= parameters.jmax / 5; j++){
            C[0][0][j] = 1;
            C[1][0][parameters.jmax - j + 1] = 1;
            C[3][parameters.jmax / 5 + j][parameters.jmax + 1] = 1;
        }
    }
}


enum obstacle_iface_orientation {
    IFACE_IS_NORTH = 1,      // upper neighbour is a fluid cell
    IFACE_IS_SOUTH = 2,      // lower neighbour is a fluid cell
    IFACE_IS_WEST  = 4,      // left  neighbour is a fluid cell
    IFACE_IS_EAST  = 8       // right neighbour is a fluid cell
};


void inner_boundary_values(
  int imax,
  int jmax,
  double **U,
  double **V,
  double **P,
  double **F,
  double **G,
  int **Flag
) {
    // loop through inner Flag field
    for (int i = 1; i <= imax; ++i) {
        for (int j = 1; j <= jmax; ++j) {
            // set no-slip conditions at obstacle cells
            switch (Flag[i][j]) {
                case IFACE_IS_NORTH:
                    // ws3 eq1.4 according to ws1 eq14-17
                    V[i  ][j]  =   0;
                    U[i-1][j]  = - U[i-1][j+1];
                    U[i  ][j]  = - U[i  ][j+1];
                    G[i  ][j]  =   V[i  ][j  ];
                    P[i  ][j]  =   P[i  ][j+1];
                    break;


                case IFACE_IS_SOUTH:
                    V[i  ][j-1]  =   0;
                    U[i-1][j  ]  = - U[i-1][j-1];
                    U[i  ][j  ]  = - U[i  ][j-1];
                    G[i  ][j-1]  =   V[i  ][j-1];
                    P[i  ][j  ]  =   P[i  ][j-1];
                    break;


                case IFACE_IS_WEST:
                    // ws3 eq1.5
                    U[i-1][j  ]  =   0;
                    V[i  ][j-1]  = - V[i-1][j-1];
                    V[i  ][j  ]  = - V[i-1][j  ];
                    F[i-1][j  ]  =   U[i-1][j  ];
                    P[i  ][j  ]  =   P[i-1][j  ];
                    break;


                case IFACE_IS_EAST:
                    U[i  ][j  ]  =   0;
                    V[i  ][j-1]  = - V[i+1][j-1];
                    V[i  ][j  ]  = - V[i+1][j  ];
                    F[i  ][j  ]  =   U[i  ][j  ];
                    P[i  ][j  ]  =   P[i+1][j  ];
                    break;


                case IFACE_IS_NORTH + IFACE_IS_EAST:
                    // IFACE_IS_NORTH:
                    V[i  ][j]  =   0;
                    U[i-1][j]  = - U[i-1][j+1];
                    //U[i  ][j]  = - U[i  ][j+1];
                    G[i  ][j]  =   V[i  ][j  ];
                    //P[i  ][j]  =   P[i  ][j+1];

                    // IFACE_IS_EAST:
                    U[i  ][j  ]  =   0;
                    V[i  ][j-1]  = - V[i+1][j-1];
                    //V[i  ][j  ]  = - V[i+1][j  ];
                    F[i  ][j  ]  =   U[i  ][j  ];
                    //P[i  ][j  ]  =   P[i+1][j  ];

                    // average pressure:
                    P[i][j] = (P[i][j+1] + P[i+1][j])/2;
                    break;


                case IFACE_IS_NORTH + IFACE_IS_WEST:
                    // IFACE_IS_NORTH:
                    V[i  ][j]  =   0;
                    //U[i-1][j]  = - U[i-1][j+1];
                    U[i  ][j]  = - U[i  ][j+1];
                    G[i  ][j]  =   V[i  ][j  ];
                    //P[i  ][j]  =   P[i  ][j+1];

                    // IFACE_IS_WEST:
                    U[i-1][j  ]  =   0;
                    V[i  ][j-1]  = - V[i-1][j-1];
                    //V[i  ][j  ]  = - V[i-1][j  ];
                    F[i-1][j  ]  =   U[i-1][j  ];
                    //P[i  ][j  ]  =   P[i-1][j  ];

                    // average pressure:
                    P[i][j] = (P[i][j+1] + P[i-1][j])/2;
                    break;


                case IFACE_IS_SOUTH + IFACE_IS_EAST:
                    // IFACE_IS_SOUTH:
                    V[i  ][j-1]  =   0;
                    U[i-1][j  ]  = - U[i-1][j-1];
                    //U[i  ][j  ]  = - U[i  ][j-1];
                    G[i  ][j-1]  =   V[i  ][j-1];
                    //P[i  ][j  ]  =   P[i  ][j-1];

                    // IFACE_IS_EAST:
                    U[i  ][j  ]  =   0;
                    //V[i  ][j-1]  = - V[i+1][j-1];
                    V[i  ][j  ]  = - V[i+1][j  ];
                    F[i  ][j  ]  =   U[i  ][j  ];
                    //P[i  ][j  ]  =   P[i+1][j  ];

                    // average pressure:
                    P[i][j] = (P[i][j-1] + P[i+1][j])/2;
                    break;


                case IFACE_IS_SOUTH + IFACE_IS_WEST:
                    // IFACE_IS_SOUTH:
                    V[i  ][j-1]  =   0;
                    //U[i-1][j  ]  = - U[i-1][j-1];
                    U[i  ][j  ]  = - U[i  ][j-1];
                    G[i  ][j-1]  =   V[i  ][j-1];
                    //P[i  ][j  ]  =   P[i  ][j-1];

                    // IFACE_IS_WEST:
                    U[i-1][j  ]  =   0;
                    //V[i  ][j-1]  = - V[i-1][j-1];
                    V[i  ][j  ]  = - V[i-1][j  ];
                    F[i-1][j  ]  =   U[i-1][j  ];
                    //P[i  ][j  ]  =   P[i-1][j  ];

                    // average pressure:
                    P[i][j] = (P[i][j-1] + P[i-1][j])/2;
                    break;
            }
        }
    }
}
