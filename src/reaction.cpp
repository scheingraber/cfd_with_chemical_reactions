#include "reaction.h"

// FInds K(T), the reaction rate temperature factos for a given pair of components
inline double compute_k(
        double activation, // Activation energy divided by R
        double freq_factor,// frequency factor
        double T // Temperature in absolute scale
        ){
    return freq_factor * exp( -activation / T );
}

// This thing calculates the reaction rate in a single point for a single
// reaction.
double single_reaction_rate(
        double ***C,            // Array of concentration matrices
        double **T,
        double T_inf,
        const reaction_t & reac,  // A reaction
        int i,
        int j
        ){

    double temp = 1;
    double value;

    // Deal with the reactants
    for ( unsigned int k = 0; k < reac.reagents.size(); k++ ){
        temp *= pow (C[reac.reagents[k]][i][j], reac.exponents_reagents[k]);
    }
    value =
        compute_k( reac.activation_E_forth, reac.freq_factor_forth, T[i][j] + T_inf )
        * temp;

    // Now with the products
    temp = 1;
    for ( unsigned int k = 0; k < reac.products.size(); k++ ){
        temp *= pow (C[reac.products[k]][i][j], reac.exponents_products[k]);
    }
    value -=
        compute_k( reac.activation_E_back, reac.freq_factor_back, T[i][j] + T_inf )
        * temp;

    return value;

}

// Fills the reaction rate vector for a single point. That is, total rates for
// every substance at a given point
void compute_reaction_rate_vector(
        double ***C,
        double **T,
        const Parameters & params,
        std::vector<double> & rates,
        int i,
        int j
        ){

    double rr;

    // Reset the reaction rate vector for a new point
    std::fill ( rates.begin(), rates.end(), 0 );

    for (unsigned int k = 0; k < params.reactions.size(); k++ ){

        // Get the reaction rate for reaction k
        rr = single_reaction_rate(C, T, params.T_inf, params.reactions[k], i, j);

        // 4th nested loop series! This one to store effects on every substance
        // Reagents decrease with forward reaction
        for( unsigned int s = 0; s < params.reactions[k].reagents.size(); s++ ){
            rates[params.reactions[k].reagents[s]]
                -= rr * params.reactions[k].st_coeff_reagents[s];
        }

        // Products increase with forward reaction
        for( unsigned int s = 0; s < params.reactions[k].products.size(); s++ ){
            rates[params.reactions[k].products[s]]
                += rr * params.reactions[k].st_coeff_products[s];
        }
    }
    // By now, the vector should contain the total reaction rates for everyone,
    // Indexed in the same way as the C matrices are.
}

double reaction_max_dt(
        double ***C,
        double **T,
        int** Flag,
        const Parameters & params,
        std::vector<double> & rates
        ){

    double max_dt = DBL_MAX;

    // Sweep the whole domain
    for (int i = 1; i <= params.imax; i++){
        for (int j = 1; j <= params.jmax; j++ ){

            // Check if the position is in the fluid
            if (Flag[i][j] & 16){

                // Force the concentration to be positive
                // The reason is too long, send a message if curious
                for ( unsigned int k = 0; k < rates.size(); k++ ){
                    C[k][i][j] = (C[k][i][j] + fabs(C[k][i][j])) / 2;
                }

                compute_reaction_rate_vector(C, T, params, rates, i, j);

                // Now see what would happen to every component
                for ( unsigned int k = 0; k < rates.size(); k++ ){

                    // Now I say that we only care if the product is being consumed
                    // I'm not totally sure about that, but sounds good
                    if ( rates[k] < 0 && (C[k][i][j] / -rates[k]) < max_dt ){
                        max_dt = C[k][i][j] / -rates[k];
                    }
                }

            }

            /*
            if ((C[0][i][j] && C[1][i][j]) || C[2][i][j]){
                std::cout << "i: " << i << ", j: " << j
                    << " <" << C[0][i][j] << ", " << C[1][i][j] << ", " << C[2][i][j] << ">" << std::endl
                    << "Rrate: [" << rates[0] << ", " << rates[1] << ", " << rates[2] << "]" << std::endl
                    << "Timestep: " << max_dt << std::endl << std::endl;
            }
            */

        }
    }
    assert(max_dt > 0);
    return max_dt;
}


// React!
void compute_reaction(
        double ***C,
        double **T,
        int **Flag,
        double dt,
        const Parameters & params,
        std::vector<double> & rates
        ){

    double heat_production;

    // For every point in the domain
    for (int i = 1; i <= params.imax; i++){
        for (int j = 1; j <= params.jmax; j++ ){

            // that isn't an obstacle
            if (Flag[i][j] & 16){

                heat_production = 0;

                // Get the reaction rate
                compute_reaction_rate_vector(C, T, params, rates, i, j);

                // And use an Euler integrator for every substance
                for ( unsigned int k = 0; k < params.substance.size(); k++ ){

                    C[k][i][j] += dt * rates[k];

                    // Compute heat production
                    heat_production -= params.substance[k].H_formation * rates[k];
                }

                // Translate heat production to temperature change
                T[i][j] += heat_production / params.vol_cp;
            }

        }
    }
}
