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

extern std::set<unsigned int> inflow_instream_dirs;
extern std::set<unsigned int> outflow_instream_dirs;

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
        border_swap_information &bsi,
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
        const border_swap_information &bsi,
        std::vector<double> &distribution_values,    
        const access_function access_function,
        const std::tuple<std::vector<unsigned int>, std::vector<unsigned int>> &y_values
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
    sim_data_tuple parallel_ts_stream_and_collide
    (
        const std::vector<start_end_it_tuple> &fluid_nodes,
        const border_swap_information &bsi,
        std::vector<double> &distribution_values,    
        const access_function access_function,
        const std::tuple<std::vector<unsigned int>, std::vector<unsigned int>> &y_values
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

    /**
     * @brief Realizes inflow and outflow by an inward stream of each border node.
     *        This method is intended for use with two-step, swap and shift algorithms.
     * 
     * @param distribution_values a vector containing the distribution values of all nodes
     * @param access_function the access function used to access the distribution values
     */
    void ghost_stream_inout
    (
        std::vector<double> &distribution_values, 
        const access_function access_function,
        const std::tuple<std::vector<unsigned int>, std::vector<unsigned int>> &y_vals
    );

    /**
     * @brief Performs a halfway bounce-back streaming update for all fluid nodes within the simulation domain.
     *        This version utilizes the ghost nodes bordering a boundary node. It is intended for use with
     *        the two-step, swap and shift algorithms.
     * 
     * @param ba see documentation of border_adjacency
     * @param distribution_values a vector containing the distribution values of all nodes
     * @param access_function the access function used to access the distribution values
     */
    void perform_boundary_update
    (
        const border_swap_information &bsi,
        std::vector<double> &distribution_values, 
        const access_function access_function
    );
}

#endif