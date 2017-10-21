#include "uvp.h"
#include "Parameters.h"
#include "helper.h"
#include <algorithm>

/**
 * Determines the value of U and G according to the formula
 *
 * @f$ F_{i,j} := u_{i,j} + \delta t \left( \frac{1}{Re} \left( \left[
    \frac{\partial^2 u}{\partial x^2} \right]_{i,j} + \left[
    \frac{\partial^2 u}{\partial y^2} \right]_{i,j} \right) - \left[
    \frac{\partial (u^2)}{\partial x} \right]_{i,j} - \left[
    \frac{\partial (uv)}{\partial y} \right]_{i,j} + g_x \right) @f$
 *
 * @f$ i=1,\ldots,imax-1, \quad j=1,\ldots,jmax @f$
 *
 * @f$ G_{i,j} := v_{i,j} + \delta t \left( \frac{1}{Re} \left(
   \left[ \frac{\partial^2 v}{\partial x^2}\right]_{i,j} + \left[ \frac{\partial^2 v}{\partial
                   y^2} \right]_{i,j} \right) - \left[ \frac{\partial
                   (uv)}{\partial x} \right]_{i,j} - \left[
                 \frac{\partial (v^2)}{\partial y} \right]_{i,j} + g_y
               \right) @f$
 *
 * @f$ i=1,\ldots,imax, \quad j=1,\ldots,jmax-1 @f$
 *
 */

// Put these signatures here to keep the header untouched

double du2dx(int i, int j, double **U, double **V, double dx, double dy, double alpha);
double duvdy(int i, int j, double **U, double **V, double dx, double dy, double alpha);
double duvdx(int i, int j, double **U, double **V, double dx, double dy, double alpha);
double dv2dy(int i, int j, double **U, double **V, double dx, double dy, double alpha);

void calculate_fg(
  const Parameters & parameters,
  double **U,
  double **V,
  double **T,
  double **F,
  double **G,
  int **Flag,
  double dt
){
    /* Compute the rest of the values */
    for (int i = 1; i <= parameters.imax; ++i){
        for (int j = 1; j <= parameters.jmax; ++j){

            if ( Flag[i][j] & 16 ){ // If we have a fluid cell
                if ( Flag[i+1][j] & 16 ){ // If the following cell in x is fluid
                    F[i][j] =
                        U[i][j]
                        + dt * (
                                1 / parameters.Re * (
                                    ( U[i+1][j] - 2 * U[i][j] + U[i-1][j] ) / (parameters.dx * parameters.dx) +
                                    ( U[i][j+1] - 2 * U[i][j] + U[i][j-1] ) / (parameters.dy * parameters.dy) )
                                - du2dx(i, j, U, V, parameters.dx, parameters.dy, parameters.alpha)
                                - duvdy(i, j, U, V, parameters.dx, parameters.dy, parameters.alpha)
                                + parameters.GX * (1 - parameters.beta / 2) * (T[i][j] + T[i+1][j])
                               );
                }


                if ( Flag[i][j+1] & 16 ){
                    G[i][j] =
                        V[i][j]
                        + dt * (
                                1 / parameters.Re * (
                                    ( V[i+1][j] - 2 * V[i][j] + V[i-1][j] ) / (parameters.dx * parameters.dx) +
                                    ( V[i][j+1] - 2 * V[i][j] + V[i][j-1] ) / (parameters.dy * parameters.dy) )
                                - dv2dy(i, j, U, V, parameters.dx, parameters.dy, parameters.alpha)
                                - duvdx(i, j, U, V, parameters.dx, parameters.dy, parameters.alpha)
                                + parameters.GY * (1 - parameters.beta / 2) * (T[i][j] + T[i][j+1])
                               );
                }
            }
        }
    }

    /* Boundary conditions for F and G */
    for (int j = 1; j <= parameters.jmax; ++j) {
        F[0][j] = U[0][j];
        F[parameters.imax][j] = U[parameters.imax][j];
    }
    for (int i = 1; i <= parameters.imax; ++i) {
        G[i][0] = V[i][0];
        G[i][parameters.jmax] = V[i][parameters.jmax];
    }
}


/**
 * This operation computes the right hand side of the pressure poisson equation.
 * The right hand side is computed according to the formula
 *
 * @f$ rs = \frac{1}{\delta t} \left( \frac{F^{(n)}_{i,j}-F^{(n)}_{i-1,j}}{\delta x} + \frac{G^{(n)}_{i,j}-G^{(n)}_{i,j-1}}{\delta y} \right)  @f$
 *
 */
void calculate_rs(
  double dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **F,
  double **G,
  double **RS,
  int **Flag
){
    for (int i = 1; i <= imax; ++i){
        for (int j = 1; j <= jmax; ++j){
            if (Flag[i][j] & 16){ // If the cell is fluid
                RS[i][j] = 1 / dt * (
                          ( F[i][j] - F[i-1][j] ) / dx
                        + ( G[i][j] - G[i][j-1] ) / dy
                    );
            }
        }
    }
}


/**
 * Determines the maximal time step size. The time step size is restricted
 * accordin to the CFL theorem. So the final time step size formula is given
 * by
 *
 * @f$ {\delta t} := \tau \, \min\left( \frac{Re}{2}\left(\frac{1}{{\delta x}^2} + \frac{1}{{\delta y}^2}\right)^{-1},  \frac{{\delta x}}{|u_{max}|},\frac{{\delta y}}{|v_{max}|} \right) @f$
 *
 */
void calculate_dt(
  const Parameters & parameters,
  double *dt,
  double **U,
  double **V,
  double **T,
  double ***C,
  int **Flag,
  std::vector<double> & rates
) {
    /* STEP 1
     * Calculate the minimum absolute values of U_ij and V_ij
     * with i = [0..imax+1], j = [0..jmax+1] i.e. including boundaries.
     */

    // initial maximum values are initialized to be the first elements of U,V
    double umax = fabs(**U),
           vmax = fabs(**V);

    // sweep matrix, searching for larger elements
    for (int i = 1; i < (parameters.imax + 2) * (parameters.jmax + 2); ++i) {
        if (fabs(*(*U+i)) > umax)   umax = fabs(*(*U+i));
        if (fabs(*(*V+i)) > vmax)   vmax = fabs(*(*V+i));
    }


    /* STEP 2
     * We have to decide between several possible values for dt.
     * Calculate those possibilities now.
     */
    std::vector<double> possible_dt;

    possible_dt.push_back(parameters.Re/2 / (1/pow(parameters.dx, 2) + 1/pow(parameters.dy, 2)));

    // dx might be divided by umax=0 and will yield a huge negative value for c2.
    // this will then be picked as minimum, which is bollocks. prevent this:
    possible_dt.push_back((umax <= 0) ? DBL_MAX : parameters.dx / fabs(umax));

    // analogously for c3:
    possible_dt.push_back((vmax <= 0) ? DBL_MAX : parameters.dy / fabs(vmax));

    // Stability condition for energy transport
    possible_dt.push_back(parameters.Re*parameters.Pr/2 / ( 1/pow(parameters.dx, 2) + 1/pow(parameters.dy, 2)));

    // Stability condition for substance difussion
    possible_dt.push_back(
            (parameters.nof_substances() <= 0)
            ? DBL_MAX
            : 1 / ( 2 * std::max_element(parameters.substance.begin(), parameters.substance.end(), [](substance_t a, substance_t b){ return a.lambda < b.lambda; })->lambda
              * ( 1/pow(parameters.dx, 2) + 1/pow(parameters.dy, 2)) ));

    // Stability condition for substance reaction
    possible_dt.push_back( reaction_max_dt(C, T, Flag, parameters, rates) );


    /* STEP 3
     * Now pick the smallest possible value for dt.
     */

    *dt = parameters.tau * *std::min_element(possible_dt.begin(), possible_dt.end());
}


/**
 * Calculates the new velocity values according to the formula
 *
 * @f$ u_{i,j}^{(n+1)}  =  F_{i,j}^{(n)} - \frac{\delta t}{\delta x} (p_{i+1,j}^{(n+1)} - p_{i,j}^{(n+1)}) @f$
 * @f$ v_{i,j}^{(n+1)}  =  G_{i,j}^{(n)} - \frac{\delta t}{\delta y} (p_{i,j+1}^{(n+1)} - p_{i,j}^{(n+1)}) @f$
 *
 * As always the index range is
 *
 * @f$ i=1,\ldots,imax-1, \quad j=1,\ldots,jmax @f$
 * @f$ i=1,\ldots,imax, \quad j=1,\ldots,jmax-1 @f$
 *
 * @image html calculate_uv.jpg
 */
void calculate_uv(
  double dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **U,
  double **V,
  double **F,
  double **G,
  double **P,
  int **Flag
) {
    for (int i = 1; i < imax; ++i) {
        for (int j = 1; j < jmax; ++j) {
            if ( Flag[i][j] & 16 ){
                if ( Flag[i+1][j] & 16 ){
                    // eq 7 for  i=1..imax-1, j=1..jmax-1
                    U[i][j] = F[i][j] - dt / dx * (P[i+1][j] - P[i][j]);
                }

                if( Flag[i][j+1] & 16 ){
                    // eq 8 for  i=1..imax-1, j=1..jmax-1
                    V[i][j] = G[i][j] - dt / dy * (P[i][j+1] - P[i][j]);
                }
            }
        }
    }

    for (int i = 1; i < imax; ++i) {
        if ( Flag[i][jmax] & 16 ){
            // eq 7 for  i=1..imax-1, j=jmax
            U[i][jmax] = F[i][jmax] - dt / dx * (P[i+1][jmax] - P[i][jmax]);
        }
    }

    for (int j = 1; j < jmax; ++j) {
        if ( Flag[imax][j] & 16 ){
            // eq 8 for  i=imax, j=1..jmax-1
            V[imax][j] = G[imax][j] - dt / dy * (P[imax][j+1] - P[imax][j]);
        }
    }
}


/* The following functions compute derivatives as shown in equations 4 and 5*/

double du2dx(int i, int j, double **U, double **V, double dx, double dy, double alpha){
  return 1 / (dx * 4) * (
    ( ( U[i][j] + U[i+1][j] ) * ( U[i][j] + U[i+1][j] ) -
      ( U[i-1][j] + U[i][j] ) * ( U[i-1][j] + U[i][j] ) )
    + alpha *
    ( fabs( U[i][j] + U[i+1][j] ) * ( U[i][j] - U[i+1][j] ) -
      fabs( U[i-1][j] + U[i][j] ) * ( U[i-1][j] - U[i][j] ) )
    );
}

double duvdy(int i, int j, double **U, double **V, double dx, double dy, double alpha){
  return 1 / (dy * 4 ) * (
    ( ( V[i][j] + V[i+1][j] ) * ( U[i][j] + U[i][j+1] ) -
      ( V[i][j-1] + V[i+1][j-1] ) * ( U[i][j-1] + U[i][j] ) )
    + alpha *
    ( fabs( V[i][j] + V[i+1][j] ) * ( U[i][j] - U[i][j+1] ) -
      fabs( V[i][j-1] + V[i+1][j-1] ) * ( U[i][j-1] - U[i][j] ) )
    );
}

double duvdx(int i, int j, double **U, double **V, double dx, double dy, double alpha){
  return 1 / ( dx * 4 ) * (
    ( ( U[i][j] + U[i][j+1] ) * ( V[i][j] + V[i+1][j] ) -
      ( U[i-1][j] + U[i-1][j+1] ) * ( V[i-1][j] + V[i][j] ) )
    + alpha *
    ( fabs( U[i][j] + U[i][j+1] ) * ( V[i][j] - V[i+1][j] ) -
      fabs( U[i-1][j] + U[i-1][j+1] ) * ( V[i-1][j] - V[i][j] ) )
    );
}

double dv2dy(int i, int j, double **U, double **V, double dx, double dy, double alpha){
  return 1 / ( dy * 4 ) * (
    ( ( V[i][j] + V[i][j+1] ) * ( V[i][j] + V[i][j+1] ) -
      ( V[i][j-1] + V[i][j] ) * ( V[i][j-1] + V[i][j] ) )
    + alpha *
    ( fabs( V[i][j] + V[i][j+1] ) * ( V[i][j] - V[i][j+1] ) -
      fabs( V[i][j-1] + V[i][j] ) * ( V[i][j-1] - V[i][j] ) )
    );
}
