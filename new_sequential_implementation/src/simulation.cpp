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
    border_swap_information swap_info,
    access_function access_function
)
{
    std::cout << "Entered setup of example domain " << std::endl;
    distribution_values.assign(TOTAL_NODE_COUNT * DIRECTION_COUNT, 0.5);
    //nodes.assign(TOTAL_NODE_COUNT, 0);
    for(auto i = 0; i < TOTAL_NODE_COUNT; ++i)
    {
        if(i % HORIZONTAL_NODES != 0 && i % HORIZONTAL_NODES != HORIZONTAL_NODES - 1) nodes.push_back(i);
    }
    phase_information.assign(TOTAL_NODE_COUNT, false);
    //fluid_nodes = {nodes.begin() + HORIZONTAL_NODES, nodes.end() - HORIZONTAL_NODES};
    for(auto it = nodes.begin() + HORIZONTAL_NODES; it < nodes.end() - HORIZONTAL_NODES; ++it)
    {
        fluid_nodes.push_back(*it);
        //std::cout << *it << std::endl;
    }
    
    // Make upper and lower ghost nodes solid
    for(auto x = 0; x < HORIZONTAL_NODES; ++x)
    {
        phase_information[access::get_node_index(x,0)] = true;
        phase_information[access::get_node_index(x,VERTICAL_NODES - 1)] = true;
    }
    std::cout << "Phase information set " << std::endl;

    // Setup inlet nodes
    std::vector<double> inlet_values{0,0,0,0,0,4.5,0,0,0};
    for(auto y = 1; y < VERTICAL_NODES - 2; ++y)
    {
        access::set_all_distribution_values(inlet_values, distribution_values, access::get_node_index(1, y), access_function);
    }
    std::cout << "All distribution values were set." << std::endl;

    swap_info = bounce_back::retrieve_border_swap_information(fluid_nodes, phase_information);

    std::cout << "Retrieved border swap information." << std::endl;
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
    std::cout << "Entered setup of example domain " << std::endl;
    distribution_values.assign(TOTAL_NODE_COUNT * DIRECTION_COUNT, 0.5);
    //nodes.assign(TOTAL_NODE_COUNT, 0);
    for(auto i = 0; i < TOTAL_NODE_COUNT; ++i)
    {
        if(i % HORIZONTAL_NODES != 0 && i % HORIZONTAL_NODES != HORIZONTAL_NODES - 1) nodes.push_back(i);
    }
    phase_information.assign(TOTAL_NODE_COUNT, false);
    //fluid_nodes = {nodes.begin() + HORIZONTAL_NODES, nodes.end() - HORIZONTAL_NODES};
    for(auto it = nodes.begin() + HORIZONTAL_NODES; it < nodes.end() - HORIZONTAL_NODES; ++it)
    {
        fluid_nodes.push_back(*it);
        //std::cout << *it << std::endl;
    }
    
    // Make upper and lower ghost nodes solid
    for(auto x = 0; x < HORIZONTAL_NODES; ++x)
    {
        phase_information[access::get_node_index(x,0)] = true;
        phase_information[access::get_node_index(x,VERTICAL_NODES - 1)] = true;
    }
    std::cout << "Phase information set " << std::endl;

    // Setup inlet nodes
    std::vector<double> inlet_values{0,0,0,0,0,4.5,0,0,0};
    for(auto y = 1; y < VERTICAL_NODES - 2; ++y)
    {
        access::set_all_distribution_values(inlet_values, distribution_values, access::get_node_index(1, y), access_function);
    }
    std::cout << "All distribution values were set." << std::endl;

    ba = bounce_back::retrieve_border_adjacencies(fluid_nodes, phase_information);

    std::cout << "Retrieved border adjacencies." << std::endl;
}