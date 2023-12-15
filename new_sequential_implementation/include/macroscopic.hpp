#ifndef MACROSCOPIC_HPP
#define MACROSCOPIC_HPP
#include "access.hpp"
#include "defines.hpp"

/**
 * @brief This namespace includes functions for calculating macroscopic observables of non-boundary nodes for the D2Q9I model.
 *        Notice that those functions assume an imcompressible fluid and are false otherwise!
 */
namespace macroscopic
{
    double density(vec_of_dist_val &distribution_functions);

    velocity flow_velocity(vec_of_dist_val &distribution_functions);

    void update_all_velocities(std::vector<double> &all_distributions, std::vector<velocity> &destination, access_function access_f);
}

#endif