#ifndef __TC_H_
#define __TC_H_

#include "Range2.h"

class Range2;
class Parameters;

void calculate_next_T (Range2 const &idx_range, double **U, double **V, double ***T, double ***T_new, int **flag, const Parameters & parameters, double dt);

void calculate_next_C (Range2 const &idx_range, double **U, double **V, double ***C,  double ***C_new, int **flag, const Parameters & parameters, double dt);

#endif
