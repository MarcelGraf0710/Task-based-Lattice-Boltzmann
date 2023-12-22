#ifndef TWO_STEP_SEQUENTIAL_HPP
#define TWO_STEP_SEQUENTIAL_HPP

#include <vector>
#include "defines.hpp"

namespace two_step_sequential
{
    /**
     * @brief Performs the streaming step for all fluid nodes within the simulation domain.
     * 
     * @param fluid_non_border_nodes A vector containing the indices of all fluid non-border nodes in the domain
     * @param ba see documentation of border_adjacency
     * @param values a vector containing all distribution values
     * @param access_function the access to node values will be performed according to this access function.
     */
    void perform_fast_stream
    (
        std::vector<unsigned int> &fluid_non_border_nodes, 
        border_adjacency &ba,
        std::vector<double> &values, 
        access_function access_function
    );

    /**
     * @brief Performs the sequential two-step algorithm for the specified number of iterations.
     * 
     * @param fluid_nodes A vector containing the indices of all fluid nodes in the domain
     * @param values the vector containing the distribution values of all nodes
     * @param access_function the access function according to which the values are to be accessed
     * @param iterations this many iterations will be performed
     */
    std::vector<sim_data_tuple> run
    (  
        std::vector<unsigned int> &fluid_nodes,       
        std::vector<double> &values, 
        border_adjacency &ba,
        access_function access_function,
        unsigned int iterations
    );
}

#endif