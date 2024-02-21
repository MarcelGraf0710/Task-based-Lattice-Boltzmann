#ifndef SHIFT_SEQUENTIAL_NEW_HPP
#define SHIFT_SEQUENTIAL_NEW_HPP

#include <vector>
#include <set>
#include "defines.hpp"
#include "access.hpp"
#include <iostream>
#include "collision.hpp"
#include "utils.hpp"
#include "boundaries.hpp"

#define SHIFT_OFFSET HORIZONTAL_NODES + 1

namespace shift_sequential
{
    /**
     * @brief Performs the steaming step in all directions for the fluid node with 
     *        the specified index.
     * 
     * @param source distribution values will be taken from this vector
     * @param destination distribution values will be rearranged in this vector
     * @param access_function function that will be used to access the distribution values
     * @param fluid_node the index of the node for which the streaming step is performed
     */
    inline void shift_stream
    (
        std::vector<double> &distribution_values, 
        const access_function &access_function, 
        const unsigned int fluid_node,
        const unsigned int read_offset,
        const unsigned int write_offset
    )
    {
        // std::cout << "Accessing shift_stream for node_index " << fluid_node << ", with read offset " << read_offset << " and write offset " << write_offset << std::endl;
        // std::cout << "Distribution values of fluid node BEFORE streaming (index is " << fluid_node << ")" << std::endl;
        // std::vector<double> dist_vals_debug = access::get_distribution_values_of(distribution_values, fluid_node, access_function);
        // to_console::print_vector(dist_vals_debug, dist_vals_debug.size());

        for (const auto direction : {0, 1, 2, 3, 4, 5, 6, 7, 8})
        {
            // std::cout << "Performing stream (idx = " << fluid_node + write_offset << ", dir = " << direction << ") [" << distribution_values[access_function(fluid_node + write_offset, direction)] << "] <- ";
            // std::cout << "(idx = " << access::get_neighbor(fluid_node + read_offset, invert_direction(direction)) << ", dir = " << direction 
            //           << ") [" << distribution_values[access_function(access::get_neighbor(fluid_node + read_offset, invert_direction(direction)), direction)] << "]" << std::endl;

            distribution_values[access_function(fluid_node + write_offset, direction)] =
                distribution_values[
                    access_function(
                        access::get_neighbor(fluid_node + read_offset, invert_direction(direction)), 
                        direction)];
        }

        // std::cout << "Distribution values of fluid node AFTER streaming (index is " << fluid_node + write_offset << ")" << std::endl;
        // dist_vals_debug = access::get_distribution_values_of(distribution_values, fluid_node + write_offset, access_function);
        // to_console::print_vector(dist_vals_debug, dist_vals_debug.size());
        // std::cout << std::endl;
    }

    /**
     * @brief Performs a combined collision and streaming step for the specified fluid node.
     * 
     * @param values the vector containing the distribution values of all nodes
     * @param access_function the access function according to which the values are to be accessed
     * @param node_index index of the fluid node whose values are to be manipulated
     * @param read_offset read offset of the current iteration (even: 0 / odd: N_e)
     * @param write_offset write offset of the current iteration (even: N_e / odd: 0)
     */

    /**
     * @brief 
     * 
     * @param distribution_values 
     * @param bsi 
     * @param access_function 
     * @param read_offset 
     * @param write_offset 
     */
    sim_data_tuple perform_shift_stream_and_collide_debug
    (
        std::vector<double> &distribution_values, 
        const std::vector<unsigned int> &fluid_nodes,
        const border_swap_information &bsi,
        const access_function access_function,
        const unsigned int iteration
    );

    /**
     * @brief 
     * 
     * @param distribution_values 
     * @param bsi 
     * @param access_function 
     * @param read_offset 
     * @param write_offset 
     */
    sim_data_tuple perform_shift_stream_and_collide
    (
        std::vector<double> &distribution_values, 
        const std::vector<unsigned int> &fluid_nodes,
        const border_swap_information &bsi,
        const access_function access_function,
        const unsigned int iteration
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
        std::vector<unsigned int> &fluid_nodes,       
        std::vector<double> &values, 
        border_swap_information &bsi,
        access_function access_function,
        unsigned int iterations
    );

    /**
     * @brief Retrieves an improved version of the border swap information data structure.
     *        This method does not consider inlet and outlet ghost nodes when performing bounce-back
     *        as the inserted values will be overwritten by inflow and outflow values anyways.
     * 
     * @param fluid_nodes 
     * @param phase_information 
     * @return border_swap_information 
     */
    border_swap_information retrieve_fast_border_swap_info
    (
        const std::vector<unsigned int> &fluid_nodes, 
        const std::vector<bool> &phase_information
    );

    /**
     * @brief Performs an outstream step for all border nodes in the directions where they border non-inout ghost nodes.
     *        The distribution values will be stored in the ghost nodes in inverted order such that
     *        after this method is executed, the border nodes can be treated like regular nodes when performing an instream.
     * 
     * @param bsi 
     * @param distribution_values 
     * @param access_function 
     */
    void emplace_bounce_back_values
    (
        const border_swap_information &bsi,
        std::vector<double> &distribution_values,
        const access_function access_function,
        const unsigned int read_offset
    );

    /**
     * @brief Updates the ghost nodes that represent inlet and outlet edges.
     *        When updating, a velocity border condition will be considered for the input
     *        and a density border condition for the output.
     *        The corresponding values are constants defined in "../include/"defines.hpp".
     * 
     * @param distribution_values a vector containing the distribution values of all nodes
     * @param access_function the access function used to access the distribution values
     */
    void update_velocity_input_density_output
    (
        std::vector<double> &distribution_values, 
        const access_function access_function,
        const unsigned int offset
    );

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
        border_swap_information &swap_info,
        const access_function access_function,
        const bool enable_debug
    );
}

#endif