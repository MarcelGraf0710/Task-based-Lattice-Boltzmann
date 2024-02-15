#include "../include/collision.hpp"
#include "../include/access.hpp"
#include "../include/utils.hpp"
#include <iostream>

/**
 * @brief Performs the collision step for a node with the specified distribution values, velocity and density.
 * 
 * @param values a vector containing the distribution values of this node.
 * @param u the flow velocity at this node
 * @param density the density at this node
 * @return a vector containing the updated distribution values.
 */
std::vector<double> collision::collide_bgk
(
    const std::vector<double> &values, 
    const velocity &u, 
    double density
)
{
    double sum = 0;
    // std::cout << "\033[31mAccessing collide_bgk\033[0m" << std::endl;
    // std::cout << "distribution values: " << std::endl;
    // to_console::print_vector(values, 10);
    // std::cout << "velocity: (" << u[0] << ", " << u[1] << ")" << std::endl;
    // std::cout << "Density: " << density;
    // std::cout << std::endl;


    std::vector<double> result = maxwell_boltzmann_distribution(u, density);
    // std::cout << "Corresponding Maxwell Boltzmann distribution is ";
    // to_console::print_vector(result, 10);

    for(auto i = 0; i < DIRECTION_COUNT; ++i)
    {
        result[i] = -(1/RELAXATION_TIME) * (values[i] - result[i]) + values[i];
    }

    // std::cout << "Final result is ";
    // to_console::print_vector(result, 10);
    // std::cout << "The density is ";
    for(auto current : result)
    {
       sum += current;
    }
    // std::cout << sum << std::endl;
    // std::cout << std::endl;
    return result;
}

/**
 * @brief Performs the collision step for all fluid nodes.
 *        This function is intended to be used in the two-step algorithm as the streaming and collision steps cannot be fused there.
 * 
 * @param fluid_nodes a vector containing the node indices of all fluid nodes
 * @param values a vector containing the distribution values of all nodes
 * @param all_velocities a vector containing the velocity values of all fluid nodes
 * @param all_densities a vector containing the density values of all fluid nodes
 * @param access this function is used to access the distribution values
 */
void collision::collide_all_bgk
(
    const std::vector<unsigned int> &fluid_nodes,
    std::vector<double> &values, 
    const std::vector<velocity> &all_velocities, 
    const std::vector<double> &all_densities,
    const access_function access
)
{
    for (const auto fluid_node : fluid_nodes)
    {
        std::vector<double> current_dist_values = access::get_distribution_values_of(values, fluid_node, access);
        std::vector<double> new_distributions = collision::collide_bgk(current_dist_values, all_velocities[fluid_node], all_densities[fluid_node]);
        access::set_distribution_values_of(new_distributions, values, fluid_node, access);
    }
}