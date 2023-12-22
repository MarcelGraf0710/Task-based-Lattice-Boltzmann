#include <vector>
#include "defines.hpp"
#ifndef SHIFT_SEQUENTIAL_HPP
#define SHIFT_SEQUENTIAL_HPP

namespace swap_sequential
{
    /**
     * @brief Sets up all distribution function values for boundary nodes
     * 
     * @param values the vector containing the distribution values of all nodes
     * @param access_function the access function according to which the values are to be accessed
     * @param offset the offset specific to the current iteration of the shift algorithm (even: 0 / odd: N_e)
     */
    void perform_shift_boundary_treatment
    (
        std::vector<double> values, 
        access_function access_function,
        unsigned int offset
    );

    /**
     * @brief Performs a combined collision and streaming step for the specified fluid node.
     * 
     * @param values the vector containing the distribution values of all nodes
     * @param access_function the access function according to which the values are to be accessed
     * @param node_index index of the fluid node whose values are to be manipulated
     * @param read_offset read offset of the current iteration (even: 0 / odd: N_e)
     * @param write_offset write offset of the current iteration (even: N_e / odd: 0)
     */
    void perform_shift_stream_and_collide
    (
        std::vector<double> values, 
        access_function access_function,
        unsigned int node_index,
        unsigned int read_offset,
        unsigned int write_offset
    );

    /**
     * @brief Performs the sequential shift algorithm for the specified number of iterations.
     * 
     * @param fluid_nodes A vector containing the indices of all fluid nodes in the domain
     * @param values the vector containing the distribution values of all nodes
     * @param access_function the access function according to which the values are to be accessed
     * @param iterations this many iterations will be performed
     */
    void run
    (  
        std::vector<unsigned int> fluid_nodes,       
        std::vector<double> values, 
        access_function access_function,
        unsigned int iterations
    );
}

#endif