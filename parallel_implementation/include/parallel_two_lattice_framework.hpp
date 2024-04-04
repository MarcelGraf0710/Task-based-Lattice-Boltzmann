#ifndef PARALLEL_TWO_LATTICE_FRAMEWORK_HPP
#define PARALLEL_TWO_LATTICE_FRAMEWORK_HPP

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
 * @brief This namespace contains all methods for the framework of the parallel two-lattice algorithm.
 *        Notice that the framework itself is the same for all algorithms but the respective executions need adaptions.
 */
namespace parallel_two_lattice_framework
{
    /**
     * @brief Performs the framework-based parallel two-lattice algorithm for the specified number of iterations.
     * 
     * @param fluid_nodes A vector of tuples of iterators pointing at the first and last fluid node of each domain
     * @param boundary_nodes see documentation of border_swap_information
     * @param distribution_values_0 source for even time steps and destination for odd time steps
     * @param distribution_values_1 source for odd time steps and destination for even time steps
     * @param access_function the access function according to which distribution values are to be accessed
     * @param iterations this many iterations will be performed
     */
    void run
    (  
        const std::vector<start_end_it_tuple> &fluid_nodes,       
        const border_swap_information &boundary_nodes,
        std::vector<double> &distribution_values_0, 
        std::vector<double> &distribution_values_1,   
        const access_function access_function,
        const unsigned int iterations
    );

    /**
     * @brief This method is a serialized debug version of perform_tl_stream_and_collide_parallel.
     *        It acts as a proof-of-concept method that is suitable for testing such that errors related to
     *        the framework itself rather than the actual parallelization can be spotted.
     * 
     * @param fluid_nodes A vector of tuples of iterators pointing at the first and last fluid node of each domain
     * @param bsi see documentation of border_swap_information
     * @param distribution_values_0 source for even time steps and destination for odd time steps
     * @param distribution_values_1 source for odd time steps and destination for even time steps
     * @param access_function the access function according to which distribution values are to be accessed
     * @return see documentation of sim_data_tuple
     */
    sim_data_tuple perform_tl_stream_and_collide_debug
    (
        const std::vector<start_end_it_tuple> &fluid_nodes,
        const border_swap_information &bsi,
        std::vector<double> &source, 
        std::vector<double> &destination,    
        const access_function access_function
    );

    /**
     * @brief Performs the combined streaming and collision step for all fluid nodes within the simulation domain.
     *        This method is adapted for the framework-based parallel two-lattice algorithm.
     *        The border conditions are enforced through ghost nodes.
     * 
     * @param fluid_nodes A vector of tuples of iterators pointing at the first and last fluid node of each domain
     * @param bsi see documentation of border_swap_information
     * @param distribution_values_0 source for even time steps and destination for odd time steps
     * @param distribution_values_1 source for odd time steps and destination for even time steps
     * @param access_function the access function according to which distribution values are to be accessed
     * @return see documentation of sim_data_tuple
     */
    sim_data_tuple perform_tl_stream_and_collide_parallel
    (
        const std::vector<start_end_it_tuple> &fluid_nodes,
        const border_swap_information &bsi,
        std::vector<double> &source, 
        std::vector<double> &destination,    
        const access_function access_function
    );

    /** 
     * @brief This helper function of parallel_two_lattice_framework::perform_tl_stream_and_collide_parallel
     *        is used in the HPX loop. It performs the actual streaming and collision.
     */
    void tl_stream_and_collide_helper
    (
        std::vector<double> &source, 
        std::vector<double> &destination, 
        const access_function &access_function, 
        std::vector<velocity> &velocities, 
        std::vector<double> &densities,
        const start_end_it_tuple bounds
    );
}

#endif