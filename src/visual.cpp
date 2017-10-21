#include "helper.h"
#include "visual.h"
#include "Range2.h"
#include <iostream>
#include <fstream>
#include <boost/algorithm/string/replace.hpp>

void write_scalars_double (std::ofstream &file, std::string name, double **m, Range2 &range)
{
    file << "SCALARS " << name << " float 1" << std::endl;
    file << "LOOKUP_TABLE default" << std::endl;
    for(int j = range.j.low; j <= range.j.high; j++) {
        for(int i = range.i.low; i <= range.i.high; i++) {
            file << m[i][j] << std::endl;
        }
    }
    file << std::endl;
}

void write_vtkHeader(std::ofstream &file, int imax, int jmax, 
                      double dx, double dy) {
    file << "# vtk DataFile Version 2.0" << std::endl
         << "TUM CFDlab SS2012 Carlos, Chris, Benedikt" << std::endl
         << "ASCII" << std::endl
         << std::endl
         << "DATASET STRUCTURED_GRID" << std::endl
         << "DIMENSIONS "  << imax+1 << " " << jmax+1 << " 1" << std::endl
         << "POINTS " << (imax+1)*(jmax+1) << " float" << std:: endl
         << std::endl;
}


void write_vtkPointCoordinates(std::ofstream &file, int imax, int jmax, double dx, double dy)
{
    double originX = 0.0;
    double originY = 0.0;

    for(int j = 0; j < jmax+1; j++) {
        for(int i = 0; i < imax+1; i++) {
            file << originX+(i*dx) << " " << originY+(j*dy) << " 0" << std::endl;
        }
    }
    file << std::endl;
}


void write_vtkFile(std::string const &problem,
                 int    timeStepNumber,
                 Parameters const &params,
                 double **U,
                 double **V,
                 double **P,
                 double **T,
                 double ***C) {

    int i,j;

    std::string filename = problem + "." + std::to_string(timeStepNumber) + ".vtk";

    std::ofstream file(filename);
    if (not file.good()) {
        ERROR(std::string("Failed to open " + filename).c_str());
        return;
    }

    Range2 range(1, params.imax, 1, params.jmax);

    write_vtkHeader(file, params.imax, params.jmax, params.dx, params.dy);
    write_vtkPointCoordinates(file, params.imax, params.jmax, params.dx, params.dy);

    file << "POINT_DATA " << std::to_string((params.imax+1)*(params.jmax+1)) << std::endl;
    file << std::endl;

    file << "VECTORS velocity float" << std::endl;
    for(j = 0; j < params.jmax+1; j++) {
        for(i = 0; i < params.imax+1; i++) {
            file << std::to_string((U[i][j] + U[i][j+1]) * 0.5) << " " << std::to_string((V[i][j] + V[i+1][j]) * 0.5) << " 0" << std::endl;
        }
    }
    file << std::endl;

    file << "CELL_DATA " << std::to_string(params.imax * params.jmax) << std::endl;

    write_scalars_double (file, "pressure", P, range);

    write_scalars_double (file, "temperature", T, range);

    for (unsigned int s = 0; s < params.nof_substances(); ++s) {
        std::string name = params.substance[s].name;
        boost::algorithm::replace_all(name, " ", "_");
        write_scalars_double (file, name, C[s], range);
    }

    file.close();
}


