#include "tc.h"
#include "Parameters.h"
#include "Range2.h"
#include "matrix.h"
#include "boundary_conditions.h"
#include <math.h>

// Functions to approximate derivatives. From Griebels' book, page 133 eq. 9.21

double uX_x(int i, int j, double **U, double **X, double dx, double gamma){
    return 1 / ( 2 * dx ) * (
        ( U[i][j] * ( X[i][j] + X[i+1][j] )
        - U[i-1][j] * ( X[i-1][j] + X[i][j] ) )
        + gamma *
        ( fabs ( U[i][j] ) * ( X[i][j] - X[i+1][j] )
        - fabs ( U[i-1][j] ) * ( X[i-1][j] - X[i][j] ) )
      );
}

double vX_y(int i, int j, double **V, double **X, double dy, double gamma){
    return 1 / ( 2 * dy ) * (
        ( V[i][j] * ( X[i][j] + X[i][j+1] )
        - V[i][j-1] * ( X[i][j-1] + X[i][j] ) )
        + gamma *
        ( fabs ( V[i][j] ) * ( X[i][j] - X[i][j+1] )
        - fabs ( V[i][j-1] ) * ( X[i][j-1] - X[i][j] ) )
      );
}

inline double X_xx(int i, int j, double **X, double dx){
    return (X[i+1][j] - 2*X[i][j] + X[i-1][j]) / (dx * dx);
}

inline double X_yy(int i, int j, double **X, double dy){
    return (X[i][j+1] - 2*X[i][j] + X[i][j-1]) / (dy * dy);
}

void calculate (Range2 const &idx_range, double **U, double **V, double **X, double **X_new, int **flag, Parameters const &parameters, double dt, double coeff, double production_coeff, int obstacle_type, double obstacle_value)
{
    for (int i = idx_range.i.low; i <= idx_range.i.high; ++i) {
        for (int j = idx_range.j.low; j <= idx_range.j.high; ++j) {
            // derived from [Gr98, 9.20]

            if (flag[i][j] & 16){ // This only for fluid cells
                // We don't want to change T, so we write to another matrix
                X_new[i][j] = X[i][j] + dt * (
                        - uX_x(i, j, U, X, parameters.dx, parameters.gamma)
                        - vX_y(i, j, V, X, parameters.dy, parameters.gamma)
                        + (X_xx(i, j, X, parameters.dx) + X_yy(i, j, X, parameters.dy)) / coeff
                        );
            }

            // Else, if boundary obstacle
            else if ( !(flag[i][j] & 16) && flag[i][j] ){

                int counter = 0; // Counts the number of fluid cells around

                // For the time being, only adiabatic boundaries if Neumann

                // Always dealing with average, similar that with pressure
                // but also working with Dirichlet boundaries

                // Reset boundaries before adding
                X_new[i][j] = 0;

                // First check in vertical direction
                for ( int l = -1; l <= 1; l++ ){
                    if ( flag[i][j+l] ){
                        if ( obstacle_type == boundary_condition["dirichlet"] ){ // If fixed temp
                            X_new[i][j] += 2*obstacle_value - X[i][j+l];
                        }
                        else if ( obstacle_type == boundary_condition["neumann"] ){ // If isolation
                            X_new[i][j] += X[i][j+l];
                        }
                        counter++;
                    }
                }

                // Then the horizontal direction
                for ( int l = -1; l <= 1; l++ ){
                    if ( flag[i+l][j] ){
                        if ( obstacle_type == boundary_condition["dirichlet"] ){ // If fixed temp
                            X_new[i][j] += 2*obstacle_value - X[i+l][j];
                        }
                        else if ( obstacle_type == boundary_condition["neumann"] ){ // If isolation
                            X_new[i][j] += X[i+l][j];
                        }
                        counter++;
                    }
                }

                // Then average the whole thing
                X_new[i][j] /= counter;

            }

            // Else, inner boundary. Set it to given temperature
            // This doesn't make sense with Neumann boundaries, but whatever
            else {
                X_new[i][j] = obstacle_value;
            }
        }
    }
}

void calculate_next_T (Range2 const &idx_range, double **U, double **V, double ***T, double ***T_new, int **flag, const Parameters & parameters, double dt)
{
    calculate (idx_range, U, V, *T, *T_new, flag, parameters, dt, parameters.Re * parameters.Pr, 1, parameters.otype, parameters.oterm);

    // Swap matrices
    swap2(T, T_new);
}

void calculate_next_C (Range2 const &idx_range, double **U, double **V, double ***C,  double ***C_new, int **flag, Parameters const &parameters, double dt)
{
    // TODO get the arguments right: coefficient
    for (unsigned int s = 0; s < parameters.nof_substances(); ++s) {
        calculate (idx_range, U, V, C[s], *C_new, flag, parameters,
                dt, 1 / (parameters.substance[s].lambda), 1,
                boundary_condition["neumann"], 0);
        swap2(&(C[s]), C_new);
    }
}
