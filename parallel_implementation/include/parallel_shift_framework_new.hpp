#ifndef PARALLEL_SHIFT_FRAMEWORK_NEW_HPP
#define PARALLEL_SHIFT_FRAMEWORK_NEW_HPP

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
#include "shift_sequential.hpp"


namespace parallel_shift_framework_new
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
        const std::vector<border_swap_information> &boundary_nodes,
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
        const std::vector<border_swap_information> &bsi,
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
    void buffer_update_odd_time_step
    (
        const std::tuple<unsigned int, unsigned int> &buffer_bounds,
        std::vector<double> &distribution_values,
        access_function access_function,
        const unsigned int buffer_offset
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
    void buffer_update_even_time_step
    (
        const std::tuple<unsigned int, unsigned int> &buffer_bounds,
        std::vector<double> &distribution_values,
        access_function access_function,
        const unsigned int buffer_offset
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
        std::vector<velocity> &velocities,
        std::vector<double> &densities, 
        const access_function access_function,
        const unsigned int offset
    );

    /**
     * @brief Prints all distribution values in to the console.
     *        They are displayed in the original order, i.e. the origin is located at the lower left corner of the printed distribution chart.
     * 
     * @param distribution_values a vector containing the distribution values of all nodes 
     * @param access_function the function used to access the distribution values
     */
    inline unsigned int determine_even_time_offset
    (
        const unsigned int node,
        const std::vector<std::tuple<unsigned int, unsigned int>> &buffer_ranges
    )
    {
        unsigned int result = 0;
        auto current = buffer_ranges.begin(); 
        //std::cout << "Current: " << std::get<0>(*current) << std::endl;
        //std::cout << "Entering " << std::endl;
        while(current < buffer_ranges.end() && node >= std::get<0>(*current))
        {
            result++;
            current++;
            //std::cout << "Current: " << std::get<0>(*current) << std::endl;
        }
        //std::cout << "Leaving" << std::endl;
        return result * (SHIFT_OFFSET);
    }

    /**
     * @brief Prints all distribution values in to the console.
     *        They are displayed in the original order, i.e. the origin is located at the lower left corner of the printed distribution chart.
     * 
     * @param distribution_values a vector containing the distribution values of all nodes 
     * @param access_function the function used to access the distribution values
     */
    inline unsigned int determine_odd_time_offset
    (
        const unsigned int node,
        const std::vector<std::tuple<unsigned int, unsigned int>> &buffer_ranges
    )
    {
        unsigned int result = 1;
        auto current = buffer_ranges.begin(); 
        //std::cout << "Current: " << std::get<0>(*current) << std::endl;
        //std::cout << "Entering " << std::endl;
        while(current < buffer_ranges.end() && node > std::get<1>(*current))
        {
            result++;
            current++;
            //std::cout << "Current: " << std::get<0>(*current) << std::endl;
        }
        //std::cout << "Leaving" << std::endl;
        return result * (SHIFT_OFFSET);
    }

    /**
     * @brief Prints all distribution values in to the console.
     *        They are displayed in the original order, i.e. the origin is located at the lower left corner of the printed distribution chart.
     * 
     * @param distribution_values a vector containing the distribution values of all nodes 
     * @param access_function the function used to access the distribution values
     */
    inline void print_distribution_values
    (
        const std::vector<double> &distribution_values, 
        const access_function access_function,
        const unsigned int offset,
        const std::vector<std::tuple<unsigned int, unsigned int>> &buffer_ranges
    )
    {
        std::vector<std::vector<unsigned int>> print_dirs = {{6,7,8}, {3,4,5}, {0,1,2}};
        unsigned int current_node_index = 0;
        unsigned int previous_direction = 0;
        std::vector<double> current_values(9,0);
        std::cout << std::setprecision(3) << std::fixed;
        unsigned int line_counter = 0;
        unsigned int natural_offset = 0;

        for(auto y = VERTICAL_NODES - 1; y >= 0; --y)
        {
            if(line_counter == SUBDOMAIN_HEIGHT)
            {
                std::cout << "\033[32m";
            }
            for(auto i = 0; i < 3; ++i)
            {
                auto current_row = print_dirs[i];
                for(auto x = 0; x < HORIZONTAL_NODES; ++x)
                {
                    //std::cout << "Currently at node with coords (" << x << ", " << y << ")" << std::endl;
                    if(x == 0 || x == 1 || x == HORIZONTAL_NODES - 1 || x == HORIZONTAL_NODES - 2) std::cout << std::setprecision(5) << std::fixed;
                    if(x == 0 && y == 0) std::cout << "\033[31m";
                    else if(x == (HORIZONTAL_NODES - 1) && y == (VERTICAL_NODES - 1)) std::cout << "\033[34m";

                    current_node_index = lbm_access::get_node_index(x, y);
                    if(offset == 0)
                        natural_offset = parallel_shift_framework_new::determine_even_time_offset(current_node_index, buffer_ranges);
                    else
                        natural_offset = parallel_shift_framework_new::determine_odd_time_offset(current_node_index, buffer_ranges);
                    //std::cout << "Node with index " << current_node_index << " received offset " << natural_offset << " and is now node " << current_node_index + natural_offset + offset << std::endl;
                    current_node_index = current_node_index + natural_offset;

                    current_values = lbm_access::get_distribution_values_of(distribution_values, current_node_index, access_function);

                    for(auto j = 0; j < 3; ++j)
                    {
                        auto direction = current_row[j];
                        std::cout << current_values[direction] << "  ";
                    }
                    std::cout << "\t";
                    if((x == 0 && y == 0) || (x == (HORIZONTAL_NODES - 1) && y == (VERTICAL_NODES - 1))) std::cout << "\033[0m";
                    if(x == 0 || x == 1 || x == HORIZONTAL_NODES - 1 || x == HORIZONTAL_NODES - 2) std::cout << std::setprecision(3) << std::fixed;
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
            std::cout << std::endl;
            if(line_counter == SUBDOMAIN_HEIGHT) line_counter = 0;
            else line_counter++;
            std::cout << "\033[0m";
        }
        std::cout << std::setprecision(3) << std::fixed;
    }

    /**
     * @brief Returns the array index of the stream layout
     * 
     * @param node the node in the simulation domain
     * @param direction the direction of the velocity vector
     * @return the index of the array storing the distribution values  
     */
    inline unsigned int stream(unsigned int node, unsigned int direction)
    {
        return (TOTAL_NODE_COUNT + BUFFER_COUNT * HORIZONTAL_NODES + SUBDOMAIN_COUNT * (SHIFT_OFFSET)) * direction + node;
    }

    /**
     * @brief Returns the array index of the bundle layout
     * 
     * @param node the node in the simulation domain
     * @param direction the direction of the velocity vector
     * @return the index of the array storing the distribution values  
     */
    inline unsigned int bundle(unsigned int node, unsigned int direction)
    {
        return 3 * (direction / 3) * (TOTAL_NODE_COUNT + BUFFER_COUNT * HORIZONTAL_NODES + SUBDOMAIN_COUNT * (SHIFT_OFFSET)) + (direction % 3) + 3 * node; 
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
    sim_data_tuple parallel_shift_stream_and_collide_actual
    (
        const std::vector<start_end_it_tuple> &fluid_nodes,
        const std::vector<border_swap_information> &bsi,
        std::vector<double> &distribution_values, 
        const access_function access_function,
        const std::tuple<std::vector<unsigned int>, std::vector<unsigned int>> &y_values,
        const std::vector<std::tuple<unsigned int, unsigned int>> &buffer_ranges,
        const unsigned int iteration
    );
}

#endif