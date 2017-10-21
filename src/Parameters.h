#ifndef PARAMETERS_WGR8XHF3
#define PARAMETERS_WGR8XHF3

#include <string>
#include <vector>
#include <boost/property_tree/ptree.hpp>


namespace pt = boost::property_tree;


struct substance_t {
    std::string name;
    double lambda;
    double H_formation;
    std::string init_file;
    double init_file_coeff;
    double init_value;
    unsigned int index;
};

struct reaction_t {
    double activation_E_back;
    double activation_E_forth;

    double freq_factor_back;
    double freq_factor_forth;

    std::vector<int> reagents; // ints since I was thinking indexes, but whatever
    std::vector<int> products;
    std::vector<int> st_coeff_reagents;
    std::vector<int> st_coeff_products;
    std::vector<double> exponents_reagents;
    std::vector<double> exponents_products;
};


class Parameters{
    public:
        int read_from_file(std::string const &filename, std::string &err_msg);
        unsigned int nof_substances() const;

    public:
        std::string problem;
        double Re;                /* reynolds number   */
        double Pr;                // Prandtl number
        double vol_cp;            // Density times heat capacity. Volumetric heat capacity
        double UI;                /* velocity x-direction */
        double VI;                /* velocity y-direction */
        double PI;                /* pressure */

        double TI;                // Temperature
        std::string TI_file;
        double TI_file_coeff;

        double beta;              // Thermal expansion coefficient
        double GX;                /* gravitation x-direction */
        double GY;                /* gravitation y-direction */
        double t_end;             /* end time */
        double xlength;           /* length of the domain x-dir.*/
        double ylength;           /* length of the domain y-dir.*/
        int imax;                 // Number of cells in x direction
        int jmax;                 // Number of cells in y direction
        double dt;                /* time step */
        double alpha;             /* uppwind differencing factor*/
        double gamma;             // same as above, temperature
        double omg;               /* relaxation factor */
        double tau;               /* safety factor for time step*/
        unsigned int  itermax;    /* max. number of iterations  */

        /* for pressure per time step */
        double eps;               /* accuracy bound for pressure*/
        std::string out_prefix;
        double out_dt;            /* time for output */

        int wlt;                // Temperature boundary type
        int wrt;
        int wtt;
        int wbt;

        double tl;                 // Temperature or flux at walls
        double tr;
        double tt;
        double tb;

        double T_inf;

        int otype;                 // Obstacle type
        double oterm;                 // Obstacle temperature

        // These are for velocity and pressure
        int wlvp;                    /* left boundary */
        int wrvp;                    /* right boundary */
        int wtvp;                    /* upper boundary */
        int wbvp;                    /* lower boundary */

        double pl;                 // Pressure value at left
        double pr;                 // Pressure value at right
        double pt;                 // Pressure value at top
        double pb;                 // Pressure value at bottom

        double dx;
        double dy;

        std::vector<substance_t> substance;

        // Reaction parameters

        std::vector <reaction_t> reactions;

        std::string geometry_file;

    private:
        void parse_params_problem    (pt::ptree const &property);
        void parse_params_time       (pt::ptree const &property);
        void parse_params_substances (pt::ptree const &property);
        void parse_params_reactions  (pt::ptree const &property);
        void parse_params_sor        (pt::ptree const &property);
        void parse_params_constants  (pt::ptree const &property);
        void parse_params_pressure   (pt::ptree const &property);
        void get_value_or_file       (pt::ptree const &tree, std::string const &what, double &value, std::string &file, double &file_coeff);
        int  elem_name_to_idx        (std::string const &name, std::vector<substance_t> const &vec);
};


#endif /* end of include guard: PARAMETERS_WGR8XHF3 */
