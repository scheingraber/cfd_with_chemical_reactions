#ifndef __RANDWERTE_H__
#define __RANDWERTE_H__

// forward decl.
class Parameters;

/**
 * The boundary values of the problem are set.
 */
void domain_boundary_values(
  const Parameters &parameters,
  double **U,
  double **V,
  double **T,
  double ***C
);

void inner_boundary_values(
  int imax,
  int jmax,
  double **U,
  double **V,
  double **P,
  double **F,
  double **G,
  int **Flag
);



/**
 * Depending on the problem, additional boundary values are applied by this function.
 */
void spec_boundary_val(
        const char *problem,
        const Parameters & parameters,
        double **U,
        double **V,
        double ***C
);


#endif
