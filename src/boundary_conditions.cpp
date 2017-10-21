#include "boundary_conditions.h"

#include <algorithm>

int parse_boundary_condition (std::string type)
{
    std::transform(type.begin(), type.end(), type.begin(), ::tolower);
    std::map<std::string,int>::iterator iter = boundary_condition.find(type);
    if (iter != boundary_condition.end())
        return iter->second;
    else
        throw "Boundary condition not found.";
}
