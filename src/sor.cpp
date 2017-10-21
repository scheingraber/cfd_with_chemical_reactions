#include "sor.h"
#include "boundary_conditions.h"
#include <math.h>

enum obstacle_iface_orientation {
    IFACE_IS_NORTH = 1,      // upper neighbour is a fluid cell
    IFACE_IS_SOUTH = 2,      // lower neighbour is a fluid cell
    IFACE_IS_WEST  = 4,      // left  neighbour is a fluid cell
    IFACE_IS_EAST  = 8       // right neighbour is a fluid cell
};

void sor(
  double omg,
  double dx,
  double dy,
  int    imax,
  int    jmax,
  double **P,
  double **RS,
  double *res,
  int    **Flag,
  int wl, int wr, int wt, int wb,// Use this to determine what kind of boundary we have
  double pl, double pr, double pt, double pb   // Pressures on the edges. Ignores if incorrect boundary type
) {
  int i,j, counter;
  double rloc;
  double coeff = omg/(2.0*(1.0/(dx*dx)+1.0/(dy*dy)));

  /* SOR iteration */
  for(i = 1; i <= imax; i++) {
    for(j = 1; j <=jmax; j++) {
      if (Flag[i][j] & 16){
        P[i][j] = (1.0-omg)*P[i][j]
                + coeff*(( P[i+1][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]+P[i][j-1])/(dy*dy) - RS[i][j]);
      }
    }
  }
  
  // Extra loop for obstacle boundaries
  for(i = 1; i <= imax; i++) {
    for(j = 1; j <=jmax; j++) {
      if ( !(Flag[i][j] & 16) && Flag[i][j]){ // If an obstacle but not an inner obstacle (Boundary cell)
          P[i][j] = 0;  // Reset pressure
          counter = 0;
          for (int l = -1; l <= 1; l+=2){
              if (Flag[i][j+l]){ // If one of the neighbours is fluid
                  P[i][j] += P[i][j+l];
                  counter++;
              }
          }
          for (int l = -1; l <= 1; l+=2){
              if (Flag[i+l][j]){ // If one of the neighbours is fluid
                  P[i][j] += P[i+l][j];
                  counter++;
              }
          }
          P[i][j] /= counter; // To get an average if it's necessary
          // The case where counter = 0 should not happen, since one of the
          // conditions is to have at least one fluid neighbour
      }
    }
  }

  /* compute the residual */
  rloc = 0;
  counter = 0;
  for(i = 1; i <= imax; i++) {
    for(j = 1; j <= jmax; j++) {
        if (Flag[i][j] & 16){
            rloc += ( (P[i+1][j]-2.0*P[i][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]-2.0*P[i][j]+P[i][j-1])/(dy*dy) - RS[i][j])*
                ( (P[i+1][j]-2.0*P[i][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]-2.0*P[i][j]+P[i][j-1])/(dy*dy) - RS[i][j]);
            counter++;
        }
    }
  }
  rloc = rloc/counter;
  rloc = sqrt(rloc);
  /* set residual */
  *res = rloc;


  // The worksheet states that we can have outflow only on the left or right,
  // so I'll just be taking those two into account.

  if ( wl == boundary_condition["pressure"] ){
      for (j = 1; j <= jmax; j++){
          P[0][j] = 2 * pl - P[1][j];
      }
  }
  else{
      for (j = 1; j <= jmax; j++){
          P[0][j] = P[1][j];
      }
  }

  if ( wr == boundary_condition["pressure"] ){
      for (j = 1; j <= jmax; j++){
          P[imax+1][j] = 2 * pr - P[imax][j];
      }
  }
  else{
      for (j = 1; j <= jmax; j++){
          P[imax+1][j] = P[imax][j];
      }
  }

  if ( wt == boundary_condition["pressure"] ){
      for (i = 1; i <= imax; i++){
          P[i][jmax+1] = 2 * pt - P[i][jmax];
      }
  }
  else{
      for (i = 1; i <= imax; i++){
          P[i][jmax+1] = P[i][jmax];
      }
  }

  if ( wb == boundary_condition["pressure"] ){
      for (i = 1; i <= imax; i++){
          P[i][0] = 2 * pt - P[i][1];
      }
  }
  else{
      for (i = 1; i <= imax; i++){
          P[i][0] = P[i][1];
      }
  }


  /* set boundary values */
  // These will always be Neumann
  for(i = 1; i <= imax; i++) {
    P[i][0] = P[i][1];
    P[i][jmax+1] = P[i][jmax];
  }

}

