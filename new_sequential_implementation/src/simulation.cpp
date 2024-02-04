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
    std::vector<double> values = {0,1,2,3,4,5,6,7,8};//maxwell_boltzmann_distribution(velocity_vectors[4], 1);
    for(auto i = 0; i < TOTAL_NODE_COUNT; ++i)
    {
        distribution_values.insert(distribution_values.end(), values.begin(), values.end());
    }
    std::cout << "Set the distribution values of all nodes to ";
    to_console::print_vector(values, DIRECTION_COUNT + 1); 
    
    // TODO: PUT THIS IN RUN!!!!
    std::vector<double> inlet_values{0,0,0,0,1,0.5,0,0,0};
    for(auto y = 1; y < VERTICAL_NODES - 1; ++y)
    {
        access::set_all_distribution_values(inlet_values, distribution_values, access::get_node_index(1, y), access_function);
    }
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
    border_adjacency &ba,
    access_function access_function
)
{
    std::cout << "Setting up example domain." << std::endl;

    /* Set up distribution values */
    std::vector<double> values = {0,0,0,0,1,0,0,0,0};
    for(auto i = 0; i < TOTAL_NODE_COUNT; ++i)
    {
        distribution_values.insert(distribution_values.end(), values.begin(), values.end());
    }
    std::cout << "Set the distribution values of all nodes to ";
    to_console::print_vector(values); 
    
    // TODO: PUT THIS IN RUN!!!!
    std::vector<double> inlet_values{0,0,0,0,1,0.5,0,0,0};
    for(auto y = 1; y < VERTICAL_NODES - 2; ++y)
    {
        access::set_all_distribution_values(inlet_values, distribution_values, access::get_node_index(1, y), access_function);
    }
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

    /* Set up boundary adjacencies */
    ba = bounce_back::retrieve_border_adjacencies(fluid_nodes, phase_information);
    std::cout << "Retrieved border adjacencies." << std::endl;
}