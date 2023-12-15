#ifndef TWO_LATTICE_SEQUENTIAL_HPP
#define TWO_LATTICE_SEQUENTIAL_HPP

#include <vector>
#include "defines.hpp"

namespace two_lattice_sequential
{
    /**
     * @brief Performs the streaming step for 
     * 
     * @param time_step current iteration of the simulation
     * @param fluid_nodes A vector containing the indices of all fluid nodes in the domain
     * @param values_0 source for even time steps and destination for odd time steps
     * @param values_1 source for odd time steps and destination for even time steps
     * @param access_function the access function according to which the values are to be accessed
     */
    void perform_tl_stream_and_collide(
        int time_step,
        std::vector<double> fluid_nodes,
        std::vector<double> values_0, 
        std::vector<double> values_1,    
        access_function access_function
    );

    /**
     * @brief Performs the sequential two-lattice algorithm for the specified number of iterations.
     * 
     * @param fluid_nodes A vector containing the indices of all fluid nodes in the domain
     * @param values the vector containing the distribution values of all nodes
     * @param access_function the access function according to which the values are to be accessed
     * @param iterations this many iterations will be performed
     */
    void run(  
        std::vector<unsigned int> fluid_nodes,       
        std::vector<double> values, 
        access_function access_function,
        unsigned int iterations
        );
}

#endif