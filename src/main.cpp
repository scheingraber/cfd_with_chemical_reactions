#include "helper.h"
#include "matrix.h"
#include "visual.h"
#include "init.h"
#include "uvp.h"
#include "pgm.h"
#include "boundary_val.h"
#include "sor.h"
#include "tc.h"
#include "reaction.h"
#include <float.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include "Range2.h"
#include "Parameters.h"

int main(int argc, char** argv){
    Parameters params;
    double dt;

    // Use a parameter file given on the command line or fallback to default
    // note: conf_dir must contain trailing slash
    std::string conf_file = "../conf/wire.xml", conf_dir  = "../conf/";
    if (argc == 2) {
        std::string arg(argv[1]);
        // c++ only way in order to avoid using boost::filesystem
        auto pivot = std::find   (arg.rbegin(), arg.rend(), '/' ).base();
        conf_file  = std::string (pivot, arg.end());
        conf_dir   = std::string (arg.begin(), pivot);
    }

    std::cout << "Reading parameters from file " << conf_file << std::endl;
    std::string err_msg;
    if (params.read_from_file(conf_dir + conf_file, err_msg) != 0) {
        ERROR(err_msg.c_str());
    }

    std::cout << "Initialising matrices..." << std::endl;

    // read from pgm file. the init_flag function also allocates storage
    // for the Flag field, since only that way we git rid of the imax and jmax
    // in the .dat file.
    int **Flag = init_flag((conf_dir + params.geometry_file).c_str(), &(params.imax), &(params.jmax) );

    Range2 idx_range(1, params.imax, 1, params.jmax);

    params.dx = params.xlength / params.imax;
    params.dy = params.ylength / params.jmax;

    // allocate storage for all matrices according to the parameters just read
    double **U = matrix<double>(0, params.imax + 1, 0, params.jmax +1);
    double **V = matrix<double>(0, params.imax + 1, 0, params.jmax +1);
    double **P = matrix<double>(0, params.imax + 1, 0, params.jmax +1);
    // Indexes to kmax and not kmax + 1 because those values are not used.
    double **F = matrix<double>(0, params.imax, 0, params.jmax);
    double **G = matrix<double>(0, params.imax, 0, params.jmax);
    double **RS = matrix<double>(0, params.imax, 0, params.jmax);

    // Concentration matrix: array of pointers to matrices of substances
    // and concentration matrix to swap
    double ***C = new double**[params.nof_substances()];

    for (unsigned int s = 0; s < params.nof_substances(); ++s) {
        // allocate and initialize a matrix for the concentration of the s'th substance
        C[s] = init (params.substance[s].init_value, conf_dir + params.substance[s].init_file, params.substance[s].init_file_coeff, params.imax, params.jmax);

        // allocate a matrix for swapping the concentration
    }

    // Swap matrix for computation of explicit quantities
    double **swap = matrix<double>(0, params.imax + 1, 0, params.jmax +1);

    // temperature
    double **T = init (params.TI, params.TI_file, params.TI_file_coeff, params.imax, params.jmax);

    // Assign initial values to u, v, p
    init_matrices(params.UI, params.VI, params.PI, params.imax, params.jmax, U, V, P);

    // A vector to accumulate the rates for each substance seems adequate.
    // Change if there's a prettier solution
    std::vector<double> rates;
    rates.resize(params.nof_substances());

    std::cout << "Starting simulation..." << std::endl;

    double t = 0;
    unsigned int n = 0;

    /* A couple of numbers to control the progress feedback of the program */
    double next_printing_time = 0;

    while (t < params.t_end) {

        if (t >= next_printing_time){
            printf("Currently at t = %f. Printing VTK\n",t);
            write_vtkFile(params.out_prefix, n, params, U, V, P, T, C);
            next_printing_time += params.out_dt;
        }

        // Set boundary values for u and v
        domain_boundary_values(params, U, V, T, C);
        inner_boundary_values(params.imax, params.jmax, U, V, P, F, G, Flag);
        spec_boundary_val(params.problem.c_str(), params, U, V, C);

        // Select dt according to (13)
        // The requirement of not changing the framework, forces us to make this decision here
        if ( params.tau > 0 ){
            calculate_dt(params, &dt, U, V, T, C, Flag, rates);
        }

        // Compute reaction effects
        compute_reaction ( C, T, Flag, dt, params, rates );

        // Compute concentration of all substances
        calculate_next_C (idx_range, U, V, C, &swap, Flag, params, dt);

        // TODO calculate reaction rate R

        // Compute temperature
        calculate_next_T (idx_range, U, V, &T, &swap, Flag, params, dt);

        // Compute F (n) and G(n) according to (9),(10),(17)
        calculate_fg(params, U, V, T, F, G, Flag, dt);

        // Compute the right-hand side rs of the pressure equation (11)
        calculate_rs(dt, params.dx, params.dy, params.imax, params.jmax, F, G, RS, Flag);

        unsigned int it = 0;
        double res = DBL_MAX;

        while ((it < params.itermax) && (res > params.eps)) {
            // Perform a SOR iteration according to (18) using the provided function and retrieve the residual res
            sor(params.omg, params.dx, params.dy, params.imax, params.jmax, P, RS, &res, Flag, params.wlvp, params.wrvp, params.wtvp, params.wbvp, params.pl, params.pr, params.pt, params.pb);
            ++it;
        }
        printf("dt: %f, current t: %f, SOR iterations: %d\n",dt,t, it);

        /* Keep the average of the pressure at zero */
        {
            double average;
            int i, j, counter=0;
            for ( i = 0 ; i <= params.imax + 1; i++){
                for ( j = 0; j <= params.jmax + 1; j++ ){
                    if( Flag[i][j] & 16 ){
                        average += P[i][j];
                        counter++;
                    }
                }
            }
            average /= counter;
            for ( i = 0 ; i <= params.imax + 1; i++){
                for ( j = 0; j <= params.jmax + 1; j++ ){
                    if( Flag[i][j] & 16 ){
                        P[i][j] -= average;
                    }
                }
            }
        }

        // Compute u(n+1) and v (n+1) according to (7),(8)
        calculate_uv(dt, params.dx, params.dy, params.imax, params.jmax, U, V, F, G, P, Flag);



        t += dt;
        ++n;
    }

    write_vtkFile(params.out_prefix, n, params, U, V, P, T, C);


    // deallocate the storage of all matrices
    free_matrix <double> (U,    0, params.imax + 1, 0, params.jmax +1);
    free_matrix <double> (V,    0, params.imax + 1, 0, params.jmax +1);
    free_matrix <double> (P,    0, params.imax + 1, 0, params.jmax +1);
    free_matrix <double> (F,    0, params.imax,     0, params.jmax);
    free_matrix <double> (G,    0, params.imax,     0, params.jmax);
    free_matrix <double> (RS,   0, params.imax,     0, params.jmax);
    free_matrix <int>    (Flag, 0, params.imax + 1, 0, params.jmax + 1);
    free_matrix <double> (swap, 0, params.imax + 1, 0, params.jmax + 1);
    free_matrix <double> (T,    0, params.imax,     0, params.jmax);

    // free concentration matrices
    // we can't simply call delete[][][], since the inner matrix was allocated
    // using malloc
    for (unsigned int i = 0; i < params.substance.size(); ++i) {
        free_matrix <double> (C[i], 0, params.imax + 1, 0, params.jmax + 1);
    }
    delete[] C;

    return 0;
}
