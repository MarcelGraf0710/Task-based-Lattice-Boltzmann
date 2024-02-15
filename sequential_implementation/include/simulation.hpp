#ifndef SIMULATION_HPP
#define SIMULATION_HPP
#include "defines.hpp"

/**
 * @brief Create an example domain for testing purposes. The domain is a rectangle with
 *        dimensions specified in the defines file where the outermost nodes are ghost nodes.
 *        The upper and lower ghost nodes are solid whereas the leftmost and rightmost columns are fluid
 *        nodes that mark the inlet and outlet respectively.
 *        Notice that all data will be written to the parameters which are assumed to be empty initially.
 * 
 * @param distribution_values a vector containing all distribution values.
 * @param nodes a vector containing all node indices, including those of solid nodes and ghost nodes.
 * @param fluid_nodes a vector containing the indices of all fluid nodes.
 * @param phase_information a vector containing the phase information of all nodes where true means solid.
 * @param swap_info see documentation of border_swap_information
 * @param access_function the domain will be prepared for access with this access function.
 */
void setup_example_domain
(
    std::vector<double> &distribution_values,
    std::vector<unsigned int> &nodes,
    std::vector<unsigned int> &fluid_nodes,
    std::vector<bool> &phase_information,
    border_swap_information &swap_info,
    access_function access_function
);


#endif