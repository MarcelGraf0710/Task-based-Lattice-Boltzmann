#ifndef PARALLEL_SHIFT_FRAMEWORK_HPP
#define PARALLEL_SHIFT_FRAMEWORK_HPP

#include <vector>
#include <set>
#include "defines.hpp"
#include "access.hpp"
#include <iostream>
#include "collision.hpp"
#include "utils.hpp"
#include "boundaries.hpp"
#include "macroscopic.hpp"
#include "parallel_framework.hpp"


namespace parallel_shift_framework
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
        std::vector<double> &distribution_values,   
        const access_function access_function,
        const unsigned int iterations
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
    sim_data_tuple perform_shift_stream_and_collide_debug
    (
        const std::vector<start_end_it_tuple> &fluid_nodes,
        const border_swap_information &bsi,
        std::vector<double> &distribution_values, 
        const access_function access_function,
        const std::tuple<std::vector<unsigned int>, std::vector<unsigned int>> &y_values,
        const std::vector<std::tuple<unsigned int, unsigned int>> &buffer_ranges,
        const unsigned int iteration
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
    sim_data_tuple parallel_shift_stream_and_collide
    (
        const std::vector<start_end_it_tuple> &fluid_nodes,
        const border_swap_information &bsi,
        std::vector<double> &distribution_values, 
        const access_function access_function,
        const std::tuple<std::vector<unsigned int>, std::vector<unsigned int>> &y_values,
        const std::vector<std::tuple<unsigned int, unsigned int>> &buffer_ranges,
        const unsigned int iteration
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
    void copy_to_buffer
    (
        const std::tuple<unsigned int, unsigned int> &buffer_bounds,
        std::vector<double> &distribution_values,
        access_function access_function,
        const unsigned int offset
    );

    inline void perform_collision
    (
        const unsigned int node,
        std::vector<double> &distribution_values, 
        const access_function &access_function, 
        std::vector<velocity> &velocities, 
        std::vector<double> &densities,
        const unsigned int write_offset
    )
    {
        std::vector<double> current_distributions = lbm_access::get_distribution_values_of(distribution_values, node + write_offset, access_function);
        velocities[node] = macroscopic::flow_velocity(current_distributions);
        densities[node] = macroscopic::density(current_distributions);
        current_distributions = collision::collide_bgk(current_distributions, velocities[node], densities[node]);
        lbm_access::set_distribution_values_of(current_distributions, distribution_values, node + write_offset, access_function);
    }
}

#endif