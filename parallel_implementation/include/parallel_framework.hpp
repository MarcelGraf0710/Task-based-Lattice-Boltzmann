#ifndef PARALLEL_FRAMEWORK_HPP
#define PARALLEL_FRAMEWORK_HPP

#include <vector>
#include <set>
#include "defines.hpp"
#include "access.hpp"
#include <iostream>
#include "collision.hpp"
#include "utils.hpp"
#include "boundaries.hpp"
#include <set>
#include <iostream>

/**
 * @brief This convenience type definition is used to declare the fluid nodes belonging to the individual subdomains.
 *        The 0th entry points towards a pointer towards the first fluid node and the 1st entry towards a pointer towards
 *        the last fluid node belonging to a certain subdomain. 
 * 
 */
typedef std::tuple<std::vector<unsigned int>::const_iterator, std::vector<unsigned int>::const_iterator> start_end_it_tuple;

/**
 * @brief This namespace contains all methods that form the basis of the frameworks used for the parallelization of
 *        the different lattice Boltzmann algorithms.
 */
namespace parallel_framework
{
    /**
     * @brief This function is used to determine the fluid nodes belonging to a certain subdomain.
     * 
     * @param subdomain the index of the subdomain (counting starts at the bottom)
     * @param fluid_nodes a vector containing all fluid nodes within the simulation domain
     * @return see documentation of start_end_it_tuple
     */
    start_end_it_tuple get_subdomain_fluid_node_pointers
    (
        const unsigned int &subdomain,
        const std::vector<unsigned int> &fluid_nodes
    );

    /**
     * @brief Returns a tuple specifying the inclusive range boundaries for the specified buffer index.
     * 
     * @return a tuple, 0th entry: start node of buffer, 1st entry: end note of buffer
     */
    std::tuple<unsigned int, unsigned int> get_buffer_node_range
    (
        const unsigned int &buffer_index
    );

    /**
     * @brief Sets up a suitable domain for parallel computation. The domain is a rectangle with
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
    void setup_parallel_domain
    (    
        std::vector<double> &distribution_values,
        std::vector<unsigned int> &nodes,
        std::vector<unsigned int> &fluid_nodes,
        std::vector<bool> &phase_information,
        access_function access_function
    );

    /**
     * @brief Retrieves a version of the border swap information data structure that is suitable for the parallel framework.
     * 
     * @param fluid_nodes a vector containing the indices of all fluid nodes within the simulation domain
     * @param phase_information a vector containing the phase information for every vector (true means solid)
     * @return border_swap_information see documentation of border_swap_information
     */
    border_swap_information retrieve_fast_border_swap_info
    (
        const std::vector<start_end_it_tuple> &fluid_node_bounds,
        const std::vector<unsigned int> &fluid_nodes,  
        const std::vector<bool> &phase_information
    );

    /**
     * @brief Performs the pre-iteration buffer initialization for the buffer with the specified boundaries.
     *        For every buffer node, the directions pointing up will be copied from the nodes below and the
     *        directions pointing down will be copied from the nodes above.
     * 
     * @param buffer_bounds a tuple containing the first and last index of the buffer
     * @param distribution_values a vector containing all distribution values
     * @param access_function this function will be used to access the distribution values
     */
    void buffer_copy_update
    (
        std::tuple<unsigned int, unsigned int> buffer_bounds,
        std::vector<double> &distribution_values,
        access_function access_function
    );
}

#endif