#include "helper.h"
#include "matrix.h"
#include "pgm.h"
#include "init.h"


/**
 * The arrays U,V and P are initialized to the constant values UI, VI and PI on
 * the whole domain including the boundary.
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
) {
    for (int i = 0; i <= (imax + 1)*(jmax + 1); ++i) {
        *(*U + i) = UI;
        *(*V + i) = VI;
        *(*P + i) = PI;
    }
}

int is_forbidden_cell(int cell) {
    if (cell & 16) { // if I am a fluid cell
        // check if there are 2 solid opposing obstacle cell neighbors
        return ((~cell & 3) == 3 || (~cell & 12) == 12);
    }
    else { // if I am a solid cell
        return ((cell & 3) == 3 || (cell & 12) == 12);
    }
}

int** init_flag(const char* pgm_file, int *imax, int *jmax) {
    int threshold = 100; // threshold grey value for pic, above we have a solid
    int **pic = NULL;
    // read pgm file
    int size[2]; // to hold the pic size returned from read_pgm
    pic = read_pgm<int>(pgm_file, size);
    #ifdef DEBUG
    printf("intiflag: Read image with size %d x %d\n",size[0], size[1]);
    #endif // DEBUG
    *imax = size[0];
    *jmax = size[1];

    // with the new design it here is finally possible to allocate storage for
    // the Flag field array
    // Flag field containing obstacle information
    int** Flag = matrix<int>(0, *imax + 1, 0, *jmax + 1);

    /* if (size[0] != imax || size[1] != jmax) { */
    /*     ERROR("PGM file dimension does not match imax, jmax dimensions from DAT file."); */
    /* } */

    // flag field: center east west south north
    // 1 means fluid, 0 means obstacle
    // e.g. 1 0 0 0 0 -> would be a fluid cell with only obstacle neighbors
    // using decimal representation -> fluid >= 16, obstacle < 16

    // Initialize the values, including boundary
    for (int i = 0; i <= (*imax)+1; ++i) {
        for (int j = 1; j <= (*jmax)+1; ++j) {
            Flag[i][j] = 0; // initialise this cell
        }
    }

    // loop through pic (w/o boundary layer) and set **Flag field boundary bits
    for (int i = 1; i <= *imax; ++i) {
        for (int j = 1; j <= *jmax; ++j) {
            if (pic[i][j] > threshold)   Flag[i][j] += 16; // cell itself is fluid
            if (pic[i+1][j] > threshold) Flag[i][j] +=  8; // east neighbor is fluid
            if (pic[i-1][j] > threshold) Flag[i][j] +=  4; // west neighbor is fluid
            if (pic[i][j-1] > threshold) Flag[i][j] +=  2; // south neighbor is fluid
            if (pic[i][j+1] > threshold) Flag[i][j] +=  1; // north neighbor is fluid
            if (is_forbidden_cell(Flag[i][j])) {
                char errmsg[256];
                sprintf(errmsg, "[%i][%i] is a forbidden cell.", i, j);
                ERROR(errmsg);
            }
        }
    }
    return Flag;
}

double **init (double const &value, std::string const &file, double const &file_coeff, const int dimx, const int dimy)
{
    double **m = 0;
    if (value < 0) {
        // no valid init_value, hence read the initial concentration from a pgm file.
        int size[2];
        m = read_pgm <double> (file.c_str(), size);
        if ((size[0] != dimx) || (size[1] != dimy)) {
            std::string err_msg = "File " +  file + " dimensions " + std::to_string(size[0]) + "x" + std::to_string(size[1]) +
                " don't match those configured " + std::to_string(dimx) + "x" + std::to_string(dimy);
            ERROR(err_msg.c_str());
        }

        // since the value range is [0..255], we have to do some scaling
        for (int i = 0; i < (dimx + 1)*(dimy+1); ++i) *(*m + i) *= file_coeff;
    }
    else {
        // initialize matrix, using init_value
        m = matrix <double> (0, dimx+1, 0, dimy+1);
        init_matrix (m, 0, dimx+1, 0, dimy+1, value);
    }

    return m;
}


