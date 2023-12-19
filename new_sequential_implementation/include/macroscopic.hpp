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
    /**
     * @brief Calculates the density of a fluid node.
     * 
     * @param distribution_functions a vector containing all distribution functions of the respective fluid node.
     */
    inline double density(const std::vector<double> &distribution_functions)
    {
        return std::accumulate(distribution_functions.begin(), distribution_functions.end(), 0);
    }

    /**
     * @brief Calculates the flow velocity of a fluid node.
     * 
     * @param distribution_functions a vector containing all distribution functions of the respective fluid node.
     * @return velocity a two-dimensional array representing the flow velocity
     */
    velocity flow_velocity(const std::vector<double> &distribution_functions);

    /**
     * @brief Calculates the velocity values for all fluid nodes in the simuation domain.
     * 
     * @param fluid_nodes A vector containing the indices of all fluid nodes within the simulation domain.
     * @param all_distributions A vector containing all distribution values. 
     * @param access_function This function is used to access the distribution values.
     */
    std::vector<velocity> calculate_all_velocities
    (
        const std::vector<unsigned int> &fluid_nodes,
        std::vector<double> &all_distributions, 
        access_function access_function
    );
}

#endif