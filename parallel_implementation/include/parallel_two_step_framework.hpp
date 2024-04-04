#ifndef PARALLEL_TWO_STEP_FRAMEWORK_HPP
#define PARALLEL_TWO_STEP_FRAMEWORK_HPP

#include <vector>
#include <set>
#include "defines.hpp"
#include "access.hpp"
#include <iostream>
#include "collision.hpp"
#include "utils.hpp"
#include "boundaries.hpp"
#include "parallel_framework.hpp"

/**
 * @brief This namespace contains all methods for the framework of the parallel two-step algorithm.
 *        Notice that the framework itself is the same for all algorithms but the respective executions need adaptions.
 */
namespace parallel_two_step_framework
{
    /**
     * @brief Performs the parallel two-step algorithm for the specified number of iterations.
     * 
     * @param fluid_nodes A vector containing the indices of all fluid nodes in the domain
     * @param distribution_values the vector containing the distribution values of all nodes
     * @param bsi see documentation of border_swap_information
     * @param access_function the access function according to which the values are to be accessed
     * @param iterations this many iterations will be performed
     */
    void run
    (  
        const std::vector<start_end_it_tuple> &fluid_nodes,       
        std::vector<double> &distribution_values, 
        const std::vector<border_swap_information> &bsis,
        access_function access_function,
        unsigned int iterations
    );

    /**
     * @brief Performs the streaming step for all fluid nodes within the simulation domain.
     *        Notice that each node is streaming outward.
     * 
     * @param fluid_nodes A vector containing the indices of all fluid nodes in the domain
     * @param distribution_values a vector containing all distribution distribution_values
     * @param access_function the access to node values will be performed according to this access function.
     */
    void perform_stream
    (
        const start_end_it_tuple fluid_node_bounds, 
        std::vector<double> &distribution_values, 
        const access_function access_function
    );

    /**
     * @brief Performs the streaming and collision step for all fluid nodes within the simulation domain.
     *        The border conditions are enforced through ghost nodes.
     * 
     * @param fluid_nodes A vector containing the indices of all fluid nodes in the domain
     * @param bsi see documentation of border_swap_information
     * @param distribution_values a vector containing all distribution values
     * @param access_function the access to node values will be performed according to this access function.
     * @return sim_data_tuple see documentation of sim_data_tuple
     */
    sim_data_tuple perform_ts_stream_and_collide
    (
        const std::vector<unsigned int> &fluid_nodes,
        const border_swap_information &bsi,
        std::vector<double> &distribution_values,    
        const access_function access_function
    );

    /**
     * @brief Performs the streaming and collision step for all fluid nodes within the simulation domain.
     *        The border conditions are enforced through ghost nodes.
     *        This variant will print several debug comments to the console.
     * 
     * @param fluid_nodes A vector containing the indices of all fluid nodes in the domain
     * @param bsi see documentation of border_swap_information
     * @param distribution_values a vector containing all distribution values
     * @param access_function the access to node values will be performed according to this access function.
     * @return sim_data_tuple see documentation of sim_data_tuple
     */
    sim_data_tuple perform_ts_stream_and_collide_debug
    (
        const std::vector<start_end_it_tuple> &fluid_nodes,
        const std::vector<border_swap_information> &bsis,
        std::vector<double> &distribution_values,    
        const access_function access_function
    );

    /** 
     * @brief This helper function of parallel_two_lattice_framework::perform_tl_stream_and_collide_parallel
     *        is used in the HPX loop. It performs the actual streaming and collision.
     */
    void perform_collision
    (
        std::vector<double> &distribution_values, 
        const access_function &access_function, 
        std::vector<velocity> &velocities, 
        std::vector<double> &densities,
        const start_end_it_tuple bounds
    );
}

#endif