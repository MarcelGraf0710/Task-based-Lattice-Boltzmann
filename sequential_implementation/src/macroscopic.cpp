#include "../include/macroscopic.hpp"
#include "../include/utils.hpp"
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
    velocity velocity_vector{0,0};
    for(int i = 0; i < DIRECTION_COUNT; ++i)
    {
        velocity_vector = VELOCITY_VECTORS.at(i);
        flow_velocity[0] += distribution_functions[i] * velocity_vector[0];
        flow_velocity[1] += distribution_functions[i] * velocity_vector[1];
    }
    return flow_velocity;
}

/**
 * @brief Calculates the velocity values for ALL nodes in the lattice.
 *        Notice that all solid nodes will automatically be assigned the velocity (0,0) regardless of their
 *        distribution values as they act as ghost nodes, and as such as containers for temporal data.
 * 
 * @param fluid_nodes A vector containing the indices of all fluid nodes within the simulation domain.
 * @param all_distributions A vector containing all distribution values. 
 * @param access_function This function is used to access the distribution values.
 */
std::vector<velocity> macroscopic::calculate_all_velocities
(
    const std::vector<unsigned int> &fluid_nodes,
    const std::vector<double> &all_distributions, 
    const access_function access_function
)
{
    std::vector<velocity> result(TOTAL_NODE_COUNT, {0,0});

    for(const auto fluid_node : fluid_nodes)
    {
        result[fluid_node] = macroscopic::flow_velocity(
                access::get_distribution_values_of(all_distributions, fluid_node, access_function));
    }
    return result;
}

/**
 * @brief Calculates the density values for ALL nodes in the lattice.
 *        Notice that for better distinction, all solid nodes will be assigned the value 
 *        std::numeric_limits<double>::max().
 * 
 * @param fluid_nodes A vector containing the indices of all fluid nodes within the simulation domain.
 * @param all_distributions A vector containing all distribution values. 
 * @param access_function This function is used to access the distribution values.
 */
std::vector<double> macroscopic::calculate_all_densities
(
    const std::vector<unsigned int> &fluid_nodes,
    const std::vector<double> &all_distributions, 
    const access_function access_function
)
{
    std::vector<double> result(TOTAL_NODE_COUNT, -1);

    for(const auto fluid_node : fluid_nodes)
    {
        result[fluid_node] = macroscopic::density(
                access::get_distribution_values_of(all_distributions, fluid_node, access_function));
    }
    return result;    
}

/**
 * @brief Returns a simulation data tuple, i.e. a tuple containing all velocities (0) and all density values (1)
 *        for all fluid nodes using the specified distribution values and access function.
 * 
 * @param fluid_nodes a vector containing the indices of all fluid vectors
 * @param all_distributions a vector containing all distribution values
 * @param access_function This function is used to access the distribution values.
 * @return a tuple containing a tuple of the velocity and density values of all fluid nodes within the 
 */
sim_data_tuple macroscopic::get_sim_data_tuple
(
    const std::vector<unsigned int> &fluid_nodes,
    const std::vector<double> &all_distributions, 
    const access_function access_function
)
{
    std::vector<velocity> velocities(fluid_nodes.size(), {0,0});
    std::vector<double> densities(fluid_nodes.size(), 0);
    
    std::vector<double> current_distibutions;
    for(const auto fluid_node : fluid_nodes)
    {
        current_distibutions = access::get_distribution_values_of(all_distributions, fluid_node, access_function);
        velocities.push_back(macroscopic::flow_velocity(current_distibutions));
        densities.push_back(macroscopic::density(current_distibutions));
    }
    return sim_data_tuple{velocities, densities};
}
