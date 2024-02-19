#include "../include/simulation.hpp"
#include "../include/access.hpp"
#include "../include/boundaries.hpp"
#include "../include/macroscopic.hpp"
#include "../include/utils.hpp"
#include <iostream>
#include <vector>
#include <numeric>

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
)
{
    std::cout << "Setting up example domain." << std::endl;

    /* Set up distribution values */

    // std::vector<double> values_0 = {0.00,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08};
    // std::vector<double> values_1 = {-0.00,0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008};

    std::vector<double> test(TOTAL_NODE_COUNT * DIRECTION_COUNT, 0);
    distribution_values = test; 

    std::vector<double> values = maxwell_boltzmann_distribution(velocity_vectors[4], 1);

    for(auto i = 0; i < TOTAL_NODE_COUNT; ++i)
    {
        // if(i % 2) access::set_distribution_values_of(values_0, distribution_values, i, access_function);
        // else access::set_distribution_values_of(values_1, distribution_values, i, access_function);

        access::set_distribution_values_of(values, distribution_values, i, access_function);
    }
    
    boundary_conditions::initialize_inout(distribution_values, access_function);

    std::cout << "All distribution values were set." << std::endl;

    /* Set all nodes for direct access */
    for(auto i = 0; i < TOTAL_NODE_COUNT; ++i)
    {
        nodes.push_back(i);
    }

    /* Set up vector containing fluid nodes within the simulation domain. */
    for(auto it = nodes.begin() + HORIZONTAL_NODES; it < nodes.end() - HORIZONTAL_NODES; ++it)
    {
        if(((*it % HORIZONTAL_NODES) != 0) && ((*it % HORIZONTAL_NODES) != (HORIZONTAL_NODES - 1))) fluid_nodes.push_back(*it);
    }
    
    /* Phase information vector */
    phase_information.assign(TOTAL_NODE_COUNT, false);
    for(auto x = 0; x < HORIZONTAL_NODES; ++x)
    {
        phase_information[access::get_node_index(x,0)] = true;
        phase_information[access::get_node_index(x,VERTICAL_NODES - 1)] = true;
    }

    /* Set up border swap information */
    swap_info = bounce_back::retrieve_border_swap_information(fluid_nodes, phase_information);
}

