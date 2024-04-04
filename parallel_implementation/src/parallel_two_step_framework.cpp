#include "../include/parallel_two_step_framework.hpp"
#include "../include/access.hpp"
#include "../include/boundaries.hpp"
#include "../include/macroscopic.hpp"
#include "../include/collision.hpp"
#include <iostream>
#include <hpx/format.hpp>
#include <hpx/future.hpp>
#include <hpx/algorithm.hpp>
#include <hpx/execution.hpp>
#include <hpx/iostream.hpp>

/**
 * @brief Performs the parallel two-step algorithm for the specified number of iterations.
 * 
 * @param fluid_nodes A vector containing the indices of all fluid nodes in the domain
 * @param distribution_values the vector containing the distribution values of all nodes
 * @param bsi see documentation of border_swap_information
 * @param access_function the access function according to which the values are to be accessed
 * @param iterations this many iterations will be performed
 */
void parallel_two_step_framework::run
(  
    const std::vector<start_end_it_tuple> &fluid_nodes,       
    std::vector<double> &distribution_values, 
    const std::vector<border_swap_information> &bsis,
    access_function access_function,
    unsigned int iterations
)
{
    to_console::print_run_greeting("Parallel two-step algorithm", iterations);

    // Initializations relevant for buffering
    std::vector<std::tuple<unsigned int, unsigned int>> buffer_ranges;
    std::vector<unsigned int> buffer_indices;
    for (auto buffer_index = 0; buffer_index < BUFFER_COUNT; ++buffer_index)
    {
        buffer_ranges.push_back(parallel_framework::get_buffer_node_range(buffer_index));
        buffer_indices.push_back(buffer_index);
    }

    std::vector<sim_data_tuple>result(
        iterations, 
        std::make_tuple(std::vector<velocity>(TOTAL_NODE_COUNT, {0,0}), std::vector<double>(TOTAL_NODE_COUNT, 0)));

    /* Parallelization framework */
    for(auto time = 0; time < iterations; ++time)
    {
        std::cout << "\033[33mIteration " << time << ":\033[0m" << std::endl;

        // Framework-based parallel two-step: combined stream and collision
        result[time] = parallel_two_step_framework::perform_ts_stream_and_collide_debug
        (fluid_nodes, bsis, distribution_values, access_function);
        std::cout << "\tFinished iteration " << time << std::endl;
    }
    to_console::buffered::print_simulation_results(result);
    std::cout << "All done, exiting simulation. " << std::endl;
}

/**
 * @brief Performs the streaming step for all fluid nodes within the simulation domain.
 *        Notice that each node is streaming outward.
 * 
 * @param fluid_nodes A vector containing the indices of all fluid nodes in the domain
 * @param distribution_values a vector containing all distribution values
 * @param access_function the access to node values will be performed according to this access function.
 */
void parallel_two_step_framework::perform_stream
(
    const start_end_it_tuple fluid_node_bounds, 
    std::vector<double> &distribution_values, 
    const access_function access_function
)
{
    /* All directions that require left-to-right and/or bottom-to-top node iteration order */
    for(auto it = std::get<0>(fluid_node_bounds); it <= std::get<1>(fluid_node_bounds); ++it)
    {
        distribution_values[access_function(lbm_access::get_neighbor(*it, 0), 0)] = distribution_values[access_function(*it, 0)];
        distribution_values[access_function(lbm_access::get_neighbor(*it, 1), 1)] = distribution_values[access_function(*it, 1)];
        distribution_values[access_function(lbm_access::get_neighbor(*it, 2), 2)] = distribution_values[access_function(*it, 2)];
        distribution_values[access_function(lbm_access::get_neighbor(*it, 3), 3)] = distribution_values[access_function(*it, 3)];
    }

    /* All directions that require right-to-left and/or top-to-bottom node iteration order */
    for(auto it = std::get<1>(fluid_node_bounds); it >= std::get<0>(fluid_node_bounds); --it)
    {
        distribution_values[access_function(lbm_access::get_neighbor(*it, 5), 5)] = distribution_values[access_function(*it, 5)];
        distribution_values[access_function(lbm_access::get_neighbor(*it, 6), 6)] = distribution_values[access_function(*it, 6)];
        distribution_values[access_function(lbm_access::get_neighbor(*it, 7), 7)] = distribution_values[access_function(*it, 7)];
        distribution_values[access_function(lbm_access::get_neighbor(*it, 8), 8)] = distribution_values[access_function(*it, 8)];
    }
}

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
sim_data_tuple parallel_two_step_framework::perform_ts_stream_and_collide_debug
(
    const std::vector<start_end_it_tuple> &fluid_nodes,
    const std::vector<border_swap_information> &bsis,
    std::vector<double> &distribution_values,    
    const access_function access_function
)
{
    std::cout << "Distribution values before stream and collide: " << std::endl;
    to_console::buffered::print_distribution_values(distribution_values, access_function);
    std::cout << std::endl;

    std::vector<velocity> velocities(TOTAL_NODE_COUNT, velocity{0,0});
    std::vector<double> densities(TOTAL_NODE_COUNT, -1);
    std::vector<double> current_distributions(DIRECTION_COUNT, 0);

    std::vector<std::tuple<unsigned int, unsigned int>> buffer_ranges;
    std::vector<unsigned int> buffer_indices;
    for (auto buffer_index = 0; buffer_index < BUFFER_COUNT; ++buffer_index)
    {
        buffer_ranges.push_back(parallel_framework::get_buffer_node_range(buffer_index));
        buffer_indices.push_back(buffer_index);
    }

    /* Streaming */
    for(auto subdomain = 0; subdomain < SUBDOMAIN_COUNT; ++subdomain)
    {
        parallel_two_step_framework::perform_stream(fluid_nodes[subdomain], distribution_values, access_function);
    }

    /* Get missing streams from buffer */
    for(auto buffer_index : buffer_indices)
    {
        parallel_framework::copy_from_buffer(buffer_ranges[buffer_index], distribution_values, access_function);
    }
    
    std::cout << "\t Distribution values after streaming:" << std::endl;
    to_console::buffered::print_distribution_values(distribution_values, access_function);
    std::cout << std::endl;

    /* Perform bounce-back using ghost nodes */
    for(auto subdomain = 0; subdomain < SUBDOMAIN_COUNT; ++subdomain)
    {
        bounce_back::perform_boundary_update(bsis[subdomain], distribution_values, access_function);
    }
    std::cout << "\t Distribution values after bounce-back update:" << std::endl;
    to_console::buffered::print_distribution_values(distribution_values, access_function);
    std::cout << std::endl;

    /* Perform inflow and outflow using ghost nodes */
    std::cout << "Performing ghost stream inout " << std::endl;
    boundary_conditions::ghost_stream_inout(distribution_values, access_function);
    std::cout << "\t Distribution values after inflow and outflow via ghost nodes:" << std::endl;
    to_console::buffered::print_distribution_values(distribution_values, access_function);
    std::cout << std::endl;

    /* Perform collision for all fluid nodes */
    for(auto subdomain = 0; subdomain < SUBDOMAIN_COUNT; ++subdomain)
    {
        parallel_two_step_framework::perform_collision(distribution_values, access_function, velocities, densities, fluid_nodes[subdomain]);
    }
    std::cout << "Distribution values after collision: " << std::endl;
    to_console::buffered::print_distribution_values(distribution_values, access_function);
    std::cout << std::endl;

    /* Update ghost nodes */
    boundary_conditions::update_velocity_input_density_output(distribution_values, velocities, densities, access_function);
    std::cout << "Distribution values after ghost node update: " << std::endl;
    to_console::buffered::print_distribution_values(distribution_values, access_function);
    std::cout << std::endl;

    // Restore buffer correctness
    for(auto buffer_index : buffer_indices)
    {
        parallel_framework::copy_to_buffer(buffer_ranges[buffer_index], distribution_values, access_function);
    }

    sim_data_tuple result{velocities, densities};

    return result;
}

// /**
//  * @brief Performs the streaming and collision step for all fluid nodes within the simulation domain.
//  *        The border conditions are enforced through ghost nodes.
//  * 
//  * @param fluid_nodes A vector containing the indices of all fluid nodes in the domain
//  * @param bsi see documentation of border_swap_information
//  * @param distribution_values a vector containing all distribution values
//  * @param access_function the access to node values will be performed according to this access function.
//  * @return sim_data_tuple see documentation of sim_data_tuple
//  */
// sim_data_tuple parallel_two_step_framework::perform_ts_stream_and_collide
// (
//     const std::vector<unsigned int> &fluid_nodes,
//     const border_swap_information &bsi,
//     std::vector<double> &distribution_values,    
//     const access_function access_function
// )
// {
//     std::vector<velocity> velocities(TOTAL_NODE_COUNT, velocity{0,0});
//     std::vector<double> densities(TOTAL_NODE_COUNT, -1);
//     std::vector<double> current_distributions(DIRECTION_COUNT, 0);

//     /* Streaming */
//     two_step_sequential::perform_fast_stream(fluid_nodes, distribution_values, access_function);

//     /* Perform bounce-back using ghost nodes */
//     bounce_back::perform_boundary_update(bsi, distribution_values, access_function);

//     /* Perform inflow and outflow using ghost nodes */
//     boundary_conditions::ghost_stream_inout(distribution_values, access_function);

//     /* Perform collision for all fluid nodes */
//     velocities = macroscopic::calculate_all_velocities(fluid_nodes, distribution_values, access_function);
//     densities = macroscopic::calculate_all_densities(fluid_nodes, distribution_values, access_function);
//     collision::collide_all_bgk(fluid_nodes, distribution_values, velocities, densities, access_function);

//     /* Update ghost nodes */
//     boundary_conditions::update_velocity_input_density_output(distribution_values, velocities, densities, access_function);

//     sim_data_tuple result{velocities, densities};

//     return result;
// }

/** 
 * @brief This helper function of parallel_two_lattice_framework::perform_tl_stream_and_collide_parallel
 *        is used in the HPX loop. It performs the actual streaming and collision.
 */
void parallel_two_step_framework::perform_collision
(
    std::vector<double> &distribution_values, 
    const access_function &access_function, 
    std::vector<velocity> &velocities, 
    std::vector<double> &densities,
    const start_end_it_tuple bounds
)
{
    std::vector<double> current_distributions(DIRECTION_COUNT, 0);
    velocity current_velocity = {0,0};
    double current_density = 0;

    for(auto it = std::get<0>(bounds); it <= std::get<1>(bounds); ++it)
    {
        // std::cout << "accessing parallel_two_step_framework::perform_collision " << std::endl;
        // std::cout << "Currently dealing with node " << *it << std::endl;
        current_distributions = 
            lbm_access::get_distribution_values_of(distribution_values, *it, access_function);

        // std::cout << "Distribution values:" << std::endl;
        // to_console::print_vector(current_distributions, 10);

        current_velocity = macroscopic::flow_velocity(current_distributions);    
        velocities[*it] = current_velocity;

        // std::cout << "Current velocity:" << std::endl;
        // std::cout << "(" << current_velocity[0] << ", " << current_velocity[1] << ")" << std::endl;

        current_density = macroscopic::density(current_distributions);
        // std::cout << "Current density: " << current_density << std::endl;
        densities[*it] = current_density;
        
        current_distributions = collision::collide_bgk(current_distributions, current_velocity, current_density);

        // std::cout << "Got distribution values " << std::endl;
        // to_console::print_vector(current_distributions, 10);
        lbm_access::set_distribution_values_of(current_distributions, distribution_values, *it, access_function);
        // std::cout << std::endl;
    }
} 