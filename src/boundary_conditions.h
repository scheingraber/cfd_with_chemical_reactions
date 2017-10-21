#ifndef BOUNDARY_CONDITIONS_H9X5XKBA
#define BOUNDARY_CONDITIONS_H9X5XKBA

#include <map>
#include <string>

static std::map<std::string,int> boundary_condition = {
    {"dirichlet", 0},
    {"neumann",   1},
    {"no-slip",   2},
    {"free-slip", 3},
    {"outflow",   4},
    {"pressure",  5},
};


int parse_boundary_condition (std::string type);

#endif /* end of include guard: BOUNDARY_CONDITIONS_H9X5XKBA */
