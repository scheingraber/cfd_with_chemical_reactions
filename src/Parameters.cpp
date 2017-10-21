#include "Parameters.h"
#include "boundary_conditions.h"
#include <boost/property_tree/xml_parser.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/optional.hpp>
#include <string>
#include <algorithm>

unsigned int Parameters::nof_substances() const
{
    return substance.size();
}

int Parameters::read_from_file(std::string const &filename, std::string &err_msg)
{
    // BEWARE
    // class members are not reset upon call of this function, calling it multiple
    // times might lead to unexpected behaviour

    pt::ptree property;

    try {
        pt::xml_parser::read_xml(filename, property);
    }
    catch (std::runtime_error &err) {
        err_msg = "Failed to parse property file: " + std::string(err.what());
        return 1;
    }
    catch (...) {
        err_msg = "Failed to parse property file " + filename;
        return 1;
    }


    try {
        parse_params_problem   (property.get_child ("problem"));
        parse_params_time      (property.get_child ("time"));
        parse_params_sor       (property.get_child ("sor"));
        parse_params_constants (property.get_child ("constants"));
        parse_params_pressure  (property.get_child ("pressure"));

        out_prefix = property.get <std::string> ("output.prefix", "data");
        out_dt     = property.get <double> ("output.dt_value");


        UI       = property.get <double> ("velocity.init.u", 0);
        VI       = property.get <double> ("velocity.init.v", 0);
        wlvp     = parse_boundary_condition(property.get <std::string> ("velocity.boundary.left.<xmlattr>.type"));
        wrvp     = parse_boundary_condition(property.get <std::string> ("velocity.boundary.right.<xmlattr>.type"));
        wtvp     = parse_boundary_condition(property.get <std::string> ("velocity.boundary.top.<xmlattr>.type"));
        wbvp     = parse_boundary_condition(property.get <std::string> ("velocity.boundary.bottom.<xmlattr>.type"));



        get_value_or_file (property, "temperature.init", TI, TI_file, TI_file_coeff);
        T_inf    = property.get <double> ("temperature.t_inf", 273.15);

        wlt      = parse_boundary_condition(property.get <std::string> ("temperature.boundary.left.<xmlattr>.type"));
        wrt      = parse_boundary_condition(property.get <std::string> ("temperature.boundary.right.<xmlattr>.type"));
        wtt      = parse_boundary_condition(property.get <std::string> ("temperature.boundary.top.<xmlattr>.type"));
        wbt      = parse_boundary_condition(property.get <std::string> ("temperature.boundary.bottom.<xmlattr>.type"));
        otype    = parse_boundary_condition(property.get <std::string> ("temperature.boundary.inner.<xmlattr>.type"));
        tl       = property.get <double> ("temperature.boundary.left");
        tr       = property.get <double> ("temperature.boundary.right");
        tt       = property.get <double> ("temperature.boundary.top");
        tb       = property.get <double> ("temperature.boundary.bottom");
        oterm    = property.get <double> ("temperature.boundary.inner");


        auto root = property.get_child_optional("substances");
        if (root) parse_params_substances (*root);

        root = property.get_child_optional("reactions");
        if (root) parse_params_reactions (*root);
    }
    catch (std::runtime_error &err) {
        err_msg = "Failed to extract values from property file: " + std::string(err.what());
        return 1;
    }
    catch (...) {
        err_msg = "Unknown error parsing property file.";
        return 1;
    }


    err_msg = "";
    return 0;
}


void Parameters::parse_params_problem (pt::ptree const &property)
{
    problem  = property.get <std::string> ("name");
    xlength  = property.get <double> ("dimensions.x");
    xlength  = property.get <double> ("dimensions.x");
    ylength  = property.get <double> ("dimensions.y");
    geometry_file = property.get <std::string> ("geometry_file");
}


void Parameters::parse_params_time (pt::ptree const &property)
{
    dt       = property.get <double> ("step");
    t_end    = property.get <double> ("max");
    tau      = property.get <double> ("tau");
}

void Parameters::parse_params_substances (pt::ptree const &property) {
    for (auto &child : property) {
        if (child.first.compare("substance")) continue;

        substance_t sub = {};

        sub.name            = child.second.get <std::string> ("name");
        sub.lambda          = child.second.get <double>      ("lambda");
        sub.H_formation     = child.second.get <double>      ("H_formation",0);

        get_value_or_file (child.second, "init", sub.init_value, sub.init_file, sub.init_file_coeff);

        substance.push_back(sub);
    }
}

int  Parameters::elem_name_to_idx (std::string const &name, std::vector<substance_t> const &vec)
{
    auto elem = std::find_if(vec.begin(), vec.end(), [name](substance_t const &elem) {return elem.name == name;});
    return (elem != vec.end()) ? std::distance(vec.begin(), elem) : -1;
}


void Parameters::parse_params_reactions (pt::ptree const &property) {
    for (auto &child : property) {
        if (child.first.compare("reaction")) continue;

        reaction_t react = {};

        react.activation_E_back  = child.second.get <double> ("activation_energy.<xmlattr>.back");
        react.activation_E_forth = child.second.get <double> ("activation_energy.<xmlattr>.forth");
        react.freq_factor_back   = child.second.get <double> ("freq_factor.<xmlattr>.back");
        react.freq_factor_forth  = child.second.get <double> ("freq_factor.<xmlattr>.forth");

        auto root_reagents = child.second.get_child("reagents");
        for (auto &reagent : root_reagents) {
            int idx = elem_name_to_idx (reagent.second.get <std::string> ("name"), substance);
            if (idx != -1) react.reagents.push_back(idx);
            else throw "Unknown reagent in reaction.";

            react.exponents_reagents.push_back (reagent.second.get <double> ("exponent"));
            react.st_coeff_reagents.push_back  (reagent.second.get <int>    ("stoichiometric_coeff"));
        }

        auto root_products = child.second.get_child("products");
        for (auto &product : root_products) {
            int idx = elem_name_to_idx (product.second.get <std::string> ("name"), substance);
            if (idx != -1) react.products.push_back(idx);
            else throw "Unknown product in reaction.";

            react.exponents_products.push_back (product.second.get <double> ("exponent"));
            react.st_coeff_products.push_back   (product.second.get <int>    ("stoichiometric_coeff"));
        }

        reactions.push_back(react);
    }
}



void Parameters::parse_params_sor (pt::ptree const &property)
{
    itermax  = property.get <double> ("itermax");
    eps      = property.get <double> ("eps");
    omg      = property.get <double> ("omega");
    alpha    = property.get <double> ("alpha");
}


void Parameters::parse_params_constants  (pt::ptree const &property)
{
    Re       = property.get <double> ("Reynolds");
    Pr       = property.get <double> ("Prandtl");
    vol_cp   = property.get <double> ("vol_cp",1);
    GX       = property.get <double> ("gravitation.x");
    GY       = property.get <double> ("gravitation.y");
    beta     = property.get <double> ("beta");
    gamma    = property.get <double> ("gamma");
}


void Parameters::parse_params_pressure (pt::ptree const &property)
{
    PI       = property.get <double> ("init", 0);
    pl       = property.get <double> ("boundary.left");
    pr       = property.get <double> ("boundary.right");
    pt       = property.get <double> ("boundary.top");
    pb       = property.get <double> ("boundary.bottom");
}


void Parameters::get_value_or_file (pt::ptree const &tree, std::string const &what, double &value, std::string &file, double &file_coeff)
{
    value = tree.get <double> (what, -1);
    if (value < 0) {
        file       = tree.get <std::string> (what + ".file");
        file_coeff = tree.get <double>      (what + ".file.<xmlattr>.coeff");
    }
}


