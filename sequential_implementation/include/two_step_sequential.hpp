#ifndef TWO_STEP_SEQUENTIAL_HPP
#define TWO_STEP_SEQUENTIAL_HPP

#include <vector>
#include "defines.hpp"
#include <set>

namespace two_step_sequential
{
    /**
     * @brief Performs the streaming step for all fluid nodes within the simulation domain.
     * 
     */
    void perform_fast_stream
    (
        const std::vector<unsigned int> &fluid_nodes, 
        std::vector<double> &values, 
        const access_function access_function
    );

    /**
     * @brief Performs the streaming and collision step for all fluid nodes within the simulation domain.
     *        The border conditions are enforced through ghost nodes.
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
     *        This variant of the streaming and collision step will print several debug comments to the console.
     */
    sim_data_tuple perform_ts_stream_and_collide_debug
    (
        const std::vector<unsigned int> &fluid_nodes,
        const border_swap_information &bsi,
        std::vector<double> &distribution_values,    
        const access_function access_function
    );

    /**
     * @brief Performs the sequential two-step algorithm for the specified number of iterations.
     * 
     * @param fluid_nodes A vector containing the indices of all fluid nodes in the domain
     * @param values the vector containing the distribution values of all nodes
     * @param access_function the access function according to which the values are to be accessed
     * @param iterations this many iterations will be performed
     */
    void run
    (  
        std::vector<unsigned int> &fluid_nodes,       
        std::vector<double> &values, 
        border_swap_information &bsi,
        access_function access_function,
        unsigned int iterations
    );

    /**
     * @brief Determines the remaining streaming option for a node based on the specified border 
     *        information vector.
     * 
     * @param current_border_info an entry of a border_swap_information object
     * @return a set containing all remaining streaming directions
     */
    std::set<unsigned int> determine_streaming_directions
    (
        const std::vector<unsigned int> &current_border_info
    );
}

#endif