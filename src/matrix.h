#ifndef MATRIX_5PM35DG0
#define MATRIX_5PM35DG0

#include "helper.h"

// implementation has to be inside the header, since the function is a template
template <typename T> T **matrix ( int nrl, int nrh, int ncl, int nch )
{
   int i;
   int nrow = nrh - nrl + 1;    /* compute number of lines */
   int ncol = nch - ncl + 1;    /* compute number of columns */

   T **pArray  = (T**) malloc((size_t)( nrow * sizeof(T*)) );
   T  *pMatrix = (T*)  malloc((size_t)( nrow * ncol * sizeof(T)));

   if( pArray  == 0)  ERROR("Storage cannot be allocated");
   if( pMatrix == 0)  ERROR("Storage cannot be allocated");

   /* first entry of the array points to the value corrected by the 
      beginning of the column */
   pArray[0] = pMatrix - ncl; 

   /* compute the remaining array entries */
   for( i = 1; i < nrow; i++ )
   {
       pArray[i] = pArray[i-1] + ncol;
   }

   /* return the value corrected by the beginning of a line */
   return pArray - nrl;
}

// implementation has to be inside the header, since the function is a template
template <typename T> void free_matrix(T **m, int nrl, int nrh, int ncl, int nch )
{
   T **pArray  = m + nrl;
   T  *pMatrix = m[nrl]+ncl;

   free( pMatrix );
   free( pArray );
}

// implementation has to be inside the header, since the function is inline
inline void swap2 (double ***a, double ***b)
{
    double **tmp = *a;
    *a = *b;
    *b = tmp;
}

// implementation has to be inside the header, since the function is inline
inline void swap3 (double ****a, double ****b)
{
    double ***tmp = *a;
    *a = *b;
    *b = tmp;
}



#endif /* end of include guard: MATRIX_5PM35DG0 */
