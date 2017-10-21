#include "Parameters.h"

#ifndef __VISUAL_H__
#define __VISUAL_H__

/**
 * Method for writing header information in vtk format. 
 * 
 * The name of the file consists of the problem name (szProblem) 
 * and of the current time step. It gets the suffix .vtk. 
 * 
 * @param szProblem      File pointer for writing info.  
 * @param timeStepNumber Number of the current time step to be printed.  
 * @param xlength Length in x-direction
 * @param ylength Length in y-direction
 * @param imax    Maximum number of entries (?) in x-direction
 * @param jmax    Maximum number of entries (?) in y-direction
 * @param dx      Mesh size in x-direction
 * @param dy      Mesh size in x-direction
 * @param U       Velocities in x-direction
 * @param V       Velocities in y-direction
 * @param P       Pressure data
 * 
 * @author Tobias Neckel
 */
void write_vtkFile(std::string const &problem,
                  int    timeStepNumber,
                  Parameters const &parameters,
                  double **U,
                  double **V,
                  double **P,
                  double **T,
                  double ***C);

#endif
