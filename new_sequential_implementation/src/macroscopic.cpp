#include "../include/macroscopic.hpp"
#include <iostream>

/**
 * @brief Calculates the flow velocity of a fluid node.
 * 
 * @param distribution_functions a vector containing all distribution functions of the respective fluid node.
 * @return velocity a two-dimensional array representing the flow velocity
 */
velocity macroscopic::flow_velocity(const std::vector<double> &distribution_functions)
{
    velocity flow_velocity{0,0};
    velocity velocity_vector;
    for(int i = 0; i < DIRECTION_COUNT; ++i)
    {
        velocity_vector = velocity_vectors[i];
        flow_velocity[0] += distribution_functions[i] * velocity_vector[0];
        flow_velocity[1] += distribution_functions[i] * velocity_vector[1];
    }
    return flow_velocity;
}

/**
 * @brief Calculates the velocity values for all fluid nodes in the simuation domain.
 * 
 * @param fluid_nodes A vector containing the indices of all fluid nodes within the simulation domain.
 * @param all_distributions A vector containing all distribution values. 
 * @param access_function This function is used to access the distribution values.
 */
std::vector<velocity> macroscopic::calculate_all_velocities
(
    const std::vector<unsigned int> &fluid_nodes,
    const std::vector<double> &all_distributions, 
    access_function access_function
)
{
    std::vector<velocity> result;
    for(auto fluid_node : fluid_nodes)
    {
        result.push_back(macroscopic::flow_velocity(access::get_all_distribution_values(all_distributions, fluid_node, access_function)));
    }
    return result;
}

