#include "../include/parallel_shift_framework_new.hpp"
#include "../include/shift_sequential.hpp"
#include "../include/access.hpp"
#include "../include/defines.hpp"
#include "../include/boundaries.hpp"
#include "../include/collision.hpp"
#include "../include/macroscopic.hpp"
#include "../include/utils.hpp"
#include <set>
#include <iostream>
#include <hpx/format.hpp>
#include <hpx/future.hpp>
#include <hpx/algorithm.hpp>
#include <hpx/execution.hpp>
#include <hpx/iostream.hpp>

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
void parallel_shift_framework_new::run
(  
    const std::vector<start_end_it_tuple> &fluid_nodes,       
    const std::vector<border_swap_information> &boundary_nodes,
    std::vector<double> &distribution_values,   
    const access_function access_function,
    const unsigned int iterations
)
{
    to_console::print_run_greeting("parallel shift algorithm", iterations);

    // Initializations relevant for buffering
    std::vector<std::tuple<unsigned int, unsigned int>> buffer_ranges;
    std::tuple<std::vector<unsigned int>, std::vector<unsigned int>> y_values;
    parallel_framework::buffer_dimension_initializations(buffer_ranges, y_values);

    std::vector<sim_data_tuple>result(
        iterations, 
        std::make_tuple(std::vector<velocity>(TOTAL_NODE_COUNT, {0,0}), std::vector<double>(TOTAL_NODE_COUNT, 0)));

    /* Parallelization framework */
    for(auto time = 0; time < iterations; ++time)
    {
        std::cout << "\033[33mIteration " << time << ":\033[0m" << std::endl;

        // Framework-based parallel two-lattice: combined stream and collision
        // result[time] = parallel_two_lattice_framework::perform_tl_stream_and_collide_parallel
        // (fluid_nodes, boundary_nodes, source, destination, access_function, y_values, buffer_ranges);
        result[time] = parallel_shift_framework_new::parallel_shift_stream_and_collide
        (fluid_nodes, boundary_nodes, distribution_values, access_function, y_values, buffer_ranges, time);
        
        std::cout << "\tFinished iteration " << time << std::endl;
    }
    to_console::buffered::print_simulation_results(result);
    std::cout << "All done, exiting simulation. " << std::endl;
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
sim_data_tuple parallel_shift_framework_new::perform_shift_stream_and_collide_debug
(
    const std::vector<start_end_it_tuple> &fluid_nodes,
    const border_swap_information &bsi,
    std::vector<double> &distribution_values, 
    const access_function access_function,
    const std::tuple<std::vector<unsigned int>, std::vector<unsigned int>> &y_values,
    const std::vector<std::tuple<unsigned int, unsigned int>> &buffer_ranges,
    const unsigned int iteration
)
{
    unsigned int current_node = 0;
    unsigned int read_offset = 0;
    unsigned int write_offset = 0;
    std::vector<double> current_distributions(DIRECTION_COUNT, 0);
    
    std::vector<velocity> velocities(TOTAL_NODE_COUNT, velocity{0,0});
    std::vector<double> densities(TOTAL_NODE_COUNT, -1);
    std::vector<double> debug_distributions((TOTAL_NODE_COUNT  + (SHIFT_OFFSET)) * DIRECTION_COUNT, 0);

    if((iteration % 2) == 0) // Even time step, dis correct
    {
        read_offset = 0;
        write_offset = (SHIFT_OFFSET);

        debug_distributions = {distribution_values.begin() + (read_offset * DIRECTION_COUNT), distribution_values.end() - (((SHIFT_OFFSET) - read_offset) * DIRECTION_COUNT)};
        std::cout << "Starting offset: " << read_offset * DIRECTION_COUNT << std::endl;
        std::cout << "Ending offset: " << (((SHIFT_OFFSET) - read_offset) * DIRECTION_COUNT) << std::endl;
        std::cout << "Distribution values before stream and collide: " << std::endl;
        to_console::buffered::print_distribution_values(distribution_values, access_function);
        std::cout << std::endl;

        // Emplace bounce-back values
        for(const auto& border_node : bsi)
        {
            shift_sequential::emplace_bounce_back_values(bsi, distribution_values, access_function, read_offset);
        }

        std::cout << "Distribution values after emplace bounce-back: " << std::endl;
        debug_distributions = {distribution_values.begin() + (read_offset * DIRECTION_COUNT), distribution_values.end() - (((SHIFT_OFFSET) - read_offset) * DIRECTION_COUNT)};
        to_console::buffered::print_distribution_values(debug_distributions, access_function);

        // Buffer update
        for(auto buffer_index = 0; buffer_index < BUFFER_COUNT; ++buffer_index)
        {
            //parallel_shift_framework_new::copy_to_buffer(buffer_ranges[buffer_index], distribution_values, access_function, read_offset);
        }

        std::cout << "Distribution values after buffer update: " << std::endl;
         debug_distributions = {distribution_values.begin() + (read_offset * DIRECTION_COUNT), distribution_values.end() - (((SHIFT_OFFSET) - read_offset) * DIRECTION_COUNT)};
        to_console::buffered::print_distribution_values(debug_distributions, access_function);

        for(auto subdomain =  SUBDOMAIN_COUNT - 1; subdomain >= 0; --subdomain)
        {
            std::cout << "Start of iteration: " << *std::get<1>(fluid_nodes[subdomain]) << std::endl;
            std::cout << "End of iteration: " << *std::get<0>(fluid_nodes[subdomain]) << std::endl;
            for(auto it = std::get<1>(fluid_nodes[subdomain]); it >= std::get<0>(fluid_nodes[subdomain]); --it)
            {
                current_node = *it;

                // Streaming
                shift_sequential::shift_stream(distribution_values, access_function, current_node, read_offset, write_offset);
            }
        }

        std::cout << "Distribution values after streaming (properly shifted): " << std::endl;
        debug_distributions = {distribution_values.begin() + write_offset * DIRECTION_COUNT, distribution_values.end() - ((SHIFT_OFFSET) - write_offset) * DIRECTION_COUNT};
        to_console::buffered::print_distribution_values(debug_distributions, access_function);

        std::cout << "Distribution values after streaming (improperly shifted): " << std::endl;
        debug_distributions = {distribution_values.begin() + read_offset * DIRECTION_COUNT, distribution_values.end() - ((SHIFT_OFFSET) - read_offset) * DIRECTION_COUNT};
        to_console::buffered::print_distribution_values(debug_distributions, access_function);

        for(auto subdomain =  SUBDOMAIN_COUNT; subdomain >= 0; --subdomain)
        {
            for(auto it = std::get<1>(fluid_nodes[subdomain]); it >= std::get<0>(fluid_nodes[subdomain]); --it)
            {
                current_node = *it;

                // Collision
                current_distributions = lbm_access::get_distribution_values_of(distribution_values, current_node + write_offset, access_function);
                velocities[current_node] = macroscopic::flow_velocity(current_distributions);
                densities[current_node] = macroscopic::density(current_distributions);
                current_distributions = collision::collide_bgk(current_distributions, velocities[current_node], densities[current_node]);
                lbm_access::set_distribution_values_of(current_distributions, distribution_values, current_node + write_offset, access_function);
            }
        }

        std::cout << "Distribution values after collision: " << std::endl;
        debug_distributions = {distribution_values.begin() + write_offset * DIRECTION_COUNT, distribution_values.end() - ((SHIFT_OFFSET) - write_offset) * DIRECTION_COUNT};
        to_console::buffered::print_distribution_values(debug_distributions, access_function);
        // to_console::buffered::print_distribution_values(distribution_values, access_function);
    }
    else
    {
        read_offset = (SHIFT_OFFSET);
        write_offset = 0;

        debug_distributions = {distribution_values.begin() + (read_offset * DIRECTION_COUNT), distribution_values.end() - (((SHIFT_OFFSET) - read_offset) * DIRECTION_COUNT)};
        std::cout << "Starting offset: " << read_offset * DIRECTION_COUNT << std::endl;
        std::cout << "Ending offset: " << (((SHIFT_OFFSET) - read_offset) * DIRECTION_COUNT) << std::endl;
        std::cout << "Distribution values before stream and collide: " << std::endl;
        to_console::buffered::print_distribution_values(debug_distributions, access_function);
        std::cout << std::endl;

        // Emplace bounce-back values
        for(const auto& border_node : bsi)
        {
            shift_sequential::emplace_bounce_back_values(bsi, distribution_values, access_function, read_offset);
        }

        std::cout << "Distribution values after streaming (properly shifted): " << std::endl;
        debug_distributions = {distribution_values.begin() + write_offset * DIRECTION_COUNT, distribution_values.end() - ((SHIFT_OFFSET) - write_offset) * DIRECTION_COUNT};
        to_console::buffered::print_distribution_values(debug_distributions, access_function);

        // Buffer update
        for(auto buffer_index = 0; buffer_index < BUFFER_COUNT; ++buffer_index)
        {
            //parallel_shift_framework_new::copy_to_buffer(buffer_ranges[buffer_index], distribution_values, access_function, read_offset);
        }

        std::cout << "Distribution values after buffer update (properly shifted): " << std::endl;
        debug_distributions = {distribution_values.begin() + write_offset * DIRECTION_COUNT, distribution_values.end() - ((SHIFT_OFFSET) - write_offset) * DIRECTION_COUNT};
        to_console::buffered::print_distribution_values(debug_distributions, access_function);

        for(auto subdomain = 0; subdomain < SUBDOMAIN_COUNT; ++subdomain)
        {
            for(auto it = std::get<0>(fluid_nodes[subdomain]); it <= std::get<1>(fluid_nodes[subdomain]); ++it)
            {
                current_node = *it;

                // Streaming
                shift_sequential::shift_stream(distribution_values, access_function, current_node, read_offset, write_offset);
            }
        }

        std::cout << "Distribution values after streaming (properly shifted): " << std::endl;
        debug_distributions = {distribution_values.begin() + write_offset * DIRECTION_COUNT, distribution_values.end() - ((SHIFT_OFFSET) - write_offset) * DIRECTION_COUNT};
        to_console::buffered::print_distribution_values(debug_distributions, access_function);

        for(auto subdomain = 0; subdomain < SUBDOMAIN_COUNT; ++subdomain)
        {
            for(auto it = std::get<0>(fluid_nodes[subdomain]); it <= std::get<1>(fluid_nodes[subdomain]); ++it)
            {
                current_node = *it;

                // Collision
                current_distributions = lbm_access::get_distribution_values_of(distribution_values, current_node + write_offset, access_function);
                velocities[current_node] = macroscopic::flow_velocity(current_distributions);
                densities[current_node] = macroscopic::density(current_distributions);
                current_distributions = collision::collide_bgk(current_distributions, velocities[current_node], densities[current_node]);
                lbm_access::set_distribution_values_of(current_distributions, distribution_values, current_node + write_offset, access_function);
            }
        }
        std::cout << "Distribution values after collision: " << std::endl;
        debug_distributions = {distribution_values.begin() + write_offset * DIRECTION_COUNT, distribution_values.end() - ((SHIFT_OFFSET) - write_offset) * DIRECTION_COUNT};
        to_console::buffered::print_distribution_values(debug_distributions, access_function);
    }

    /* Update ghost nodes OFFSET ISSUE?*/
    shift_sequential::update_velocity_input_density_output(distribution_values, access_function, write_offset);

    std::cout << "Distribution values after inout update: " << std::endl;
    debug_distributions = {distribution_values.begin() + write_offset * DIRECTION_COUNT, distribution_values.end() - ((SHIFT_OFFSET) - write_offset) * DIRECTION_COUNT};
    to_console::buffered::print_distribution_values(debug_distributions, access_function);
    
    unsigned int update_node = 0;
    for(auto x = 0; x < HORIZONTAL_NODES; x = x + HORIZONTAL_NODES - 1)
    {
        for(auto y = 1; y < VERTICAL_NODES - 1; ++y)
        {
            update_node = lbm_access::get_node_index(x,y);
            current_distributions = 
                lbm_access::get_distribution_values_of(distribution_values, update_node + write_offset, access_function);
                
            velocities[update_node] = macroscopic::flow_velocity(current_distributions);
            densities[update_node] = macroscopic::density(current_distributions);
        }
    }

        unsigned int x = 0;
        velocity v = INLET_VELOCITY;
        // std::cout << "x = " << x << std::endl;
        int y = 0;
        // update_node = lbm_access::get_node_index(x,y);
        // current_distributions = maxwell_boltzmann_distribution(v, 1);
        // lbm_access::set_distribution_values_of(current_distributions, distribution_values, update_node + read_offset, access_function);
        // std::cout << "Set value of node " << update_node + write_offset << " to equilibrium distribution " << std::endl;
        
        // y = VERTICAL_NODES - 1;
        // update_node = lbm_access::get_node_index(x,y);
        // current_distributions = maxwell_boltzmann_distribution(v, 1);
        // lbm_access::set_distribution_values_of(current_distributions, distribution_values, update_node + write_offset, access_function);
        // std::cout << "Set value of node " << update_node + write_offset << " to equilibrium distribution " << std::endl;

        x = HORIZONTAL_NODES - 1;
        v = OUTLET_VELOCITY;
        std::cout << "x = " << x << std::endl;
        y = 0;
        update_node = lbm_access::get_node_index(x,y);
        current_distributions = maxwell_boltzmann_distribution(v, 1);
        lbm_access::set_distribution_values_of(current_distributions, distribution_values, update_node + write_offset, access_function);
        std::cout << "Set value of node " << update_node + write_offset << " to equilibrium distribution " << std::endl;
        
        y = VERTICAL_NODES - 1;
        update_node = lbm_access::get_node_index(x,y);
        current_distributions = maxwell_boltzmann_distribution(v, 1);
        lbm_access::set_distribution_values_of(current_distributions, distribution_values, update_node + write_offset, access_function);
        std::cout << "Set value of node " << update_node + write_offset << " to equilibrium distribution " << std::endl;

    sim_data_tuple result{velocities, densities};
    return result;
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
sim_data_tuple parallel_shift_framework_new::parallel_shift_stream_and_collide
(
    const std::vector<start_end_it_tuple> &fluid_nodes,
    const std::vector<border_swap_information> &bsi,
    std::vector<double> &distribution_values, 
    const access_function access_function,
    const std::tuple<std::vector<unsigned int>, std::vector<unsigned int>> &y_values,
    const std::vector<std::tuple<unsigned int, unsigned int>> &buffer_ranges,
    const unsigned int iteration
)
{
    unsigned int current_node = 0;
    unsigned int read_offset = 0;
    unsigned int write_offset = 0;
    std::vector<velocity> velocities(TOTAL_NODE_COUNT, velocity{0,0});
    std::vector<double> densities(TOTAL_NODE_COUNT, -1);

    if((iteration % 2) == 0)
    {
        read_offset = 0;
        write_offset = (SHIFT_OFFSET);

        std::cout << "Beginning with iteration " << iteration << std::endl;
        std::cout << "Is even?: " << true << std::endl;
        std::cout << "Read offset: " << read_offset << std::endl;
        std::cout << "Write offset: " << write_offset << std::endl;

        std::cout << "Distributions at the start:" << std::endl;
        parallel_shift_framework_new::print_distribution_values(distribution_values, access_function, read_offset, buffer_ranges);

        // Emplace bounce-back values
        for(auto subdomain = 0; subdomain < SUBDOMAIN_COUNT; ++subdomain)
        {
            unsigned int subdomain_offset = subdomain * (SHIFT_OFFSET);
            // std::cout << "Currently working on subdomain: " << subdomain << std::endl;
            // std::cout << "Subdomain offset: " << subdomain_offset << std::endl;
            bounce_back::emplace_bounce_back_values_parallel(bsi[subdomain], distribution_values, access_function, subdomain_offset + read_offset);
        }

        std::cout << "Distribution values after bounce-back update " << std::endl;
        parallel_shift_framework_new::print_distribution_values(distribution_values, access_function, read_offset, buffer_ranges);

        // Buffer update
        for(auto buffer = 0; buffer < BUFFER_COUNT; ++buffer)
        {
            int buffer_offset = (buffer+1) * (SHIFT_OFFSET);
            parallel_shift_framework_new::buffer_update_even_time_step(buffer_ranges[buffer], distribution_values, access_function, buffer_offset);
        }

        std::cout << "Distribution values after buffer update" << std::endl;
        parallel_shift_framework_new::print_distribution_values(distribution_values, access_function, read_offset, buffer_ranges);

        for(auto subdomain = 0; subdomain < SUBDOMAIN_COUNT; ++subdomain)
        {
            unsigned int subdomain_offset = subdomain * (SHIFT_OFFSET);
            for(auto it = std::get<1>(fluid_nodes[subdomain]); it >= std::get<0>(fluid_nodes[subdomain]); --it)
            {
                shift_sequential::shift_stream(distribution_values, access_function, *it, read_offset + subdomain_offset, write_offset + subdomain_offset);
            }
        }

        std::cout << "Distribution values after streaming " << std::endl;
        parallel_shift_framework_new::print_distribution_values(distribution_values, access_function, write_offset, buffer_ranges);

        for(auto subdomain = 0; subdomain < SUBDOMAIN_COUNT; ++subdomain)
        {
            unsigned int subdomain_offset = subdomain * (SHIFT_OFFSET);
            for(auto it = std::get<1>(fluid_nodes[subdomain]); it >= std::get<0>(fluid_nodes[subdomain]); --it)
            {
                parallel_shift_framework_new::perform_collision(*it, distribution_values, access_function, velocities, densities, write_offset + subdomain_offset);
            }
        }

        std::cout << "Distribution values after collision " << std::endl;
        parallel_shift_framework_new::print_distribution_values(distribution_values, access_function, write_offset, buffer_ranges);
    }
    else
    {
        read_offset = (SHIFT_OFFSET);
        write_offset = 0;

        std::cout << "Beginning with iteration " << iteration << std::endl;
        std::cout << "Is even?: " << false << std::endl;
        std::cout << "Read offset: " << read_offset << std::endl;
        std::cout << "Write offset: " << write_offset << std::endl;

        parallel_shift_framework_new::print_distribution_values(distribution_values, access_function, read_offset, buffer_ranges);

        // Emplace bounce-back values
        for(auto subdomain = 0; subdomain < SUBDOMAIN_COUNT; ++subdomain)
        {
            unsigned int subdomain_offset = subdomain * (SHIFT_OFFSET);
            bounce_back::emplace_bounce_back_values_parallel(bsi[subdomain], distribution_values, access_function, subdomain_offset + read_offset);
        }

        std::cout << "Distribution values after bounce-back update " << std::endl;
        parallel_shift_framework_new::print_distribution_values(distribution_values, access_function, read_offset, buffer_ranges);

        // Buffer update
        for(auto buffer = 0; buffer < BUFFER_COUNT; ++buffer)
        {
            int buffer_offset = (buffer+1) * (SHIFT_OFFSET);
            parallel_shift_framework_new::buffer_update_odd_time_step(buffer_ranges[buffer], distribution_values, access_function, buffer_offset);
        }

        std::cout << "Distribution values after buffer update" << std::endl;
        parallel_shift_framework_new::print_distribution_values(distribution_values, access_function, read_offset, buffer_ranges);

        for(auto subdomain = 0; subdomain < SUBDOMAIN_COUNT; ++subdomain)
        {
            unsigned int subdomain_offset = subdomain * (SHIFT_OFFSET);
            for(auto it = std::get<0>(fluid_nodes[subdomain]); it <= std::get<1>(fluid_nodes[subdomain]); ++it)
            {
                shift_sequential::shift_stream(distribution_values, access_function, *it, read_offset + subdomain_offset, write_offset + subdomain_offset);
            }
        }

        std::cout << "Distribution values after streaming" << std::endl;
        parallel_shift_framework_new::print_distribution_values(distribution_values, access_function, write_offset, buffer_ranges);

        for(auto subdomain = 0; subdomain < SUBDOMAIN_COUNT; ++subdomain)
        {
            unsigned int subdomain_offset = subdomain * (SHIFT_OFFSET);
            for(auto it = std::get<0>(fluid_nodes[subdomain]); it <= std::get<1>(fluid_nodes[subdomain]); ++it)
            {
                parallel_shift_framework_new::perform_collision(*it, distribution_values, access_function, velocities, densities, write_offset + subdomain_offset);
            }
        }

        std::cout << "Distribution values after collision" << std::endl;
        parallel_shift_framework_new::print_distribution_values(distribution_values, access_function, write_offset, buffer_ranges);
    }

    /* Update ghost nodes OFFSET ISSUE?*/
    parallel_shift_framework_new::update_velocity_input_density_output(distribution_values, velocities, densities, access_function, write_offset);
    std::cout << "Distribution values after input update" << std::endl;
    parallel_shift_framework_new::print_distribution_values(distribution_values, access_function, write_offset, buffer_ranges);
    
    sim_data_tuple result{velocities, densities};
    return result;
}



/**
 * @brief Performs the pre-iteration buffer initialization for the buffer with the specified boundaries.
 *        For every buffer node, the directions pointing up will be copied from the nodes below and the
 *        directions pointing down will be copied from the nodes above.
 * 
 * @param buffer_bounds a tuple containing the first and last index of the buffer
 * @param distribution_values a vector containing all distribution values
 * @param access_function this function will be used to access the distribution values
 */
void parallel_shift_framework_new::buffer_update_even_time_step
(
    const std::tuple<unsigned int, unsigned int> &buffer_bounds,
    std::vector<double> &distribution_values,
    access_function access_function,
    const unsigned int buffer_offset
)
{
    unsigned int start = std::get<0>(buffer_bounds);
    unsigned int end = std::get<1>(buffer_bounds);

    for(auto buffer_node = start; buffer_node <= end; ++buffer_node)
    {
        for(auto direction : {6,7,8})
        {
            distribution_values[access_function(buffer_node + buffer_offset, direction)] = 
            distribution_values[access_function(lbm_access::get_neighbor(buffer_node, 1) + buffer_offset - (SHIFT_OFFSET), direction)];
        }
        for(auto direction : {0,1,2})
        {
            distribution_values[access_function(lbm_access::get_neighbor(buffer_node, 1) + buffer_offset - 1, direction)] = 
            distribution_values[access_function(lbm_access::get_neighbor(buffer_node, 7) + buffer_offset, direction)];
            // std::cout << "Writing distribution value " << 
            // distribution_values[access_function(lbm_access::get_neighbor(buffer_node, 7) + buffer_offset, direction)] << 
            // " to dir " << direction << " of node " << lbm_access::get_neighbor(buffer_node, 1) + buffer_offset - (SHIFT_OFFSET) << std::endl;
        }
    }
}

/**
 * @brief Performs the pre-iteration buffer initialization for the buffer with the specified boundaries.
 *        For every buffer node, the directions pointing up will be copied from the nodes below and the
 *        directions pointing down will be copied from the nodes above.
 * 
 * @param buffer_bounds a tuple containing the first and last index of the buffer
 * @param distribution_values a vector containing all distribution values
 * @param access_function this function will be used to access the distribution values
 */
void parallel_shift_framework_new::buffer_update_odd_time_step
(
    const std::tuple<unsigned int, unsigned int> &buffer_bounds,
    std::vector<double> &distribution_values,
    access_function access_function,
    const unsigned int buffer_offset
)
{
    unsigned int start = std::get<0>(buffer_bounds);
    unsigned int end = std::get<1>(buffer_bounds);

    for(auto buffer_node = start; buffer_node <= end; ++buffer_node)
    {
        for(auto direction : {6,7,8})
        {
            distribution_values[access_function(buffer_node + buffer_offset + SHIFT_OFFSET, direction)] = 
            distribution_values[access_function(lbm_access::get_neighbor(buffer_node, 1) + buffer_offset, direction)];
        }
        for(auto direction : {0,1,2})
        {
            distribution_values[access_function(buffer_node + buffer_offset, direction)] =
            distribution_values[access_function(lbm_access::get_neighbor(buffer_node, 7) + buffer_offset + SHIFT_OFFSET, direction)];
            
        }
    }
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
void parallel_shift_framework_new::setup_parallel_domain
(    
    std::vector<double> &distribution_values,
    std::vector<unsigned int> &nodes,
    std::vector<unsigned int> &fluid_nodes,
    std::vector<bool> &phase_information,
    access_function access_function
)
{
    distribution_values.assign((TOTAL_NODE_COUNT + BUFFER_COUNT * HORIZONTAL_NODES + SUBDOMAIN_COUNT * (SHIFT_OFFSET)) * DIRECTION_COUNT, 0); 
    std::vector<double> regular_values = maxwell_boltzmann_distribution(VELOCITY_VECTORS.at(4), 1);
    std::vector<double> inlet_values = maxwell_boltzmann_distribution(INLET_VELOCITY, INLET_DENSITY);
    std::vector<double> outlet_values = maxwell_boltzmann_distribution(OUTLET_VELOCITY, OUTLET_DENSITY);

    // Set up vector of all nodes for direct access
    for(auto i = 0; i < TOTAL_NODE_COUNT; ++i)
    {
        nodes.push_back(i);
    }

    for(auto subdomain = 0; subdomain < SUBDOMAIN_COUNT; ++subdomain)
    {
        unsigned int subdomain_offset = subdomain * (SHIFT_OFFSET); // subdomain * ((SHIFT_OFFSET) + HORIZONTAL_NODES);
        int last_node = 0;
        for(auto y = subdomain * SUBDOMAIN_HEIGHT + subdomain; y < (subdomain + 1) * SUBDOMAIN_HEIGHT + subdomain; ++y)
        {
            for(auto x = 0; x < HORIZONTAL_NODES; ++x)
            {
                last_node = lbm_access::get_node_index(x,y);
                if(x == 0)
                {
                    lbm_access::set_distribution_values_of(inlet_values, distribution_values, last_node + subdomain_offset, access_function);
                }
                else if(x == HORIZONTAL_NODES - 1)
                {
                    lbm_access::set_distribution_values_of(outlet_values, distribution_values, last_node + subdomain_offset, access_function);
                }
                else
                {
                    lbm_access::set_distribution_values_of(regular_values, distribution_values, last_node + subdomain_offset, access_function);
                }
            }
        }
    }

    /* Phase information vector */
    phase_information.assign(TOTAL_NODE_COUNT, false);
    for(auto x = 0; x < HORIZONTAL_NODES; ++x)
    {
        phase_information[lbm_access::get_node_index(x,0)] = true;
        phase_information[lbm_access::get_node_index(x,VERTICAL_NODES - 1)] = true;
    }

    /* Fluid nodes vector */
    for(auto y = 1; y < VERTICAL_NODES - 1; ++y)
    {
        for(auto x = 1; x < HORIZONTAL_NODES - 1; ++x)
        {
            fluid_nodes.push_back(lbm_access::get_node_index(x,y));
        }
    }
}

/**
 * @brief Updates the ghost nodes that represent inlet and outlet edges.
 *        When updating, a velocity border condition will be considered for the input
 *        and a density border condition for the output.
 *        The corresponding values are constants defined in "../include/"defines.hpp".
 * 
 * @param distribution_values a vector containing the distribution values of all nodes
 * @param access_function the access function used to access the distribution values
 */
void parallel_shift_framework_new::update_velocity_input_density_output
(
    std::vector<double> &distribution_values,
    std::vector<velocity> &velocities,
    std::vector<double> &densities, 
    const access_function access_function,
    const unsigned int offset
)
{
    for(auto subdomain = 0; subdomain < SUBDOMAIN_COUNT; ++subdomain)
    {
        unsigned int subdomain_offset = subdomain * (SHIFT_OFFSET); // subdomain * ((SHIFT_OFFSET) + HORIZONTAL_NODES);
        int last_node = 0;
        std::vector<double> current_dist_vals(DIRECTION_COUNT, 0);
        int current_border_node = 0;
        velocity v = {0,0};

        for(auto y = subdomain * SUBDOMAIN_HEIGHT + subdomain; y < (subdomain + 1) * SUBDOMAIN_HEIGHT + subdomain; ++y)
        {
            current_border_node = lbm_access::get_node_index(0,y);
            current_dist_vals = maxwell_boltzmann_distribution(INLET_VELOCITY, INLET_DENSITY);
            lbm_access::set_distribution_values_of
            (
                current_dist_vals,
                distribution_values,
                current_border_node + offset + subdomain_offset,
                access_function
            );
            velocities[current_border_node] = macroscopic::flow_velocity(current_dist_vals);
            densities[current_border_node] = macroscopic::density(current_dist_vals);

            current_border_node = lbm_access::get_node_index(HORIZONTAL_NODES - 1,y);
            v = macroscopic::flow_velocity(lbm_access::get_distribution_values_of(distribution_values, lbm_access::get_neighbor(current_border_node + offset + subdomain_offset, 3), access_function));
            current_dist_vals = maxwell_boltzmann_distribution(v, OUTLET_DENSITY);
            lbm_access::set_distribution_values_of
            (
                current_dist_vals,
                distribution_values,
                current_border_node + offset + subdomain_offset,
                access_function
            );

            velocities[current_border_node] = v;
            densities[current_border_node] = OUTLET_DENSITY;
        }
    }

        unsigned int update_node = 0;
        std::vector<double> current_distributions = maxwell_boltzmann_distribution(OUTLET_VELOCITY, OUTLET_DENSITY);
        unsigned int x =  HORIZONTAL_NODES - 1;

        int y = 0;
        update_node = lbm_access::get_node_index(x,y);
        lbm_access::set_distribution_values_of(current_distributions, distribution_values, update_node + offset, access_function);
        
        y = VERTICAL_NODES - 1;
        std::cout << "Last node I'm dealing with is original " << lbm_access::get_node_index(x,y) << std::endl;
        std::cout << "Additional term is " << (BUFFER_COUNT) * (SHIFT_OFFSET) << std::endl;
        std::cout << "BUFFER_COUNT: " << BUFFER_COUNT << std::endl;
        std::cout << "SHIFT_OFFSET: " << SHIFT_OFFSET << std::endl;
        update_node = lbm_access::get_node_index(x,y) + (BUFFER_COUNT) * (SHIFT_OFFSET);
         std::cout << "I am being transferred to index " << update_node << ", is that correct?" << std::endl;
        lbm_access::set_distribution_values_of(current_distributions, distribution_values, update_node + offset, access_function);
}