#ifndef __INIT_H_
#define __INIT_H_

/**
 * The arrays U,V and P are initialized to the constant values UI, VI and PI on
 * the whole domain.
 */
void init_matrices(
  double UI,
  double VI,
  double PI,
  int imax,
  int jmax,
  double **U,
  double **V,
  double **P
);

/**
 * Returns 1 if flag field cell is forbidden cell, 0 otherwise.
 */
int is_forbidden_cell(int Field);

/**
 * Initialise the **Flag field with the obstacle and fluid cell flags.
 * Returns int** pointer to Flag array.
 */
int** init_flag(
  const char* pgm_file,
  int *imax,
  int *jmax);

double **init (double const &value, std::string const &file, double const &file_coeff, const int dimx, const int dimy);

#endif

