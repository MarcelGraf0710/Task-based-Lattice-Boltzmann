#include "../include/parallel_swap_framework.hpp"
#include "../include/swap_sequential.hpp"
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

std::set<unsigned int> inflow_instream_dirs{2,5,8};
std::set<unsigned int> outflow_instream_dirs{0,3,6};

/**
 * @brief Performs the parallel two-step algorithm for the specified number of iterations.
 * 
 * @param fluid_nodes a vector containing the first and last element of an iterator over all fluid nodes within each subdomain
 * @param distribution_values the vector containing the distribution values of all nodes
 * @param bsi see documentation of border_swap_information
 * @param access_function the access function according to which the values are to be accessed
 * @param iterations this many iterations will be performed
 */
void parallel_swap_framework::run
(  
    const std::vector<start_end_it_tuple> &fluid_nodes,       
    std::vector<double> &distribution_values, 
    const border_swap_information &bsi,
    const access_function access_function,
    const unsigned int iterations
)
{
    to_console::print_run_greeting("parallel swap algorithm", iterations);

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

        // Framework-based parallel two-step: combined stream and collision
        result[time] = parallel_swap_framework::parallel_swap_stream_and_collide
        (fluid_nodes, bsi, distribution_values, access_function, y_values, buffer_ranges);
        std::cout << "\tFinished iteration " << time << std::endl;
    }
    to_console::buffered::print_simulation_results(result);
    std::cout << "All done, exiting simulation. " << std::endl;
}

/**
 * @brief Performs the streaming and collision step for all fluid nodes within the simulation domain.
 *        The border conditions are enforced through ghost nodes.
 *        This variant will print several debug comments to the console.
 * 
 * @param fluid_nodes a vector containing the first and last element of an iterator over all fluid nodes within each subdomain
 * @param bsi see documentation of border_swap_information
 * @param distribution_values a vector containing all distribution values
 * @param access_function the access to node values will be performed according to this access function
 * @param y_values a tuple containing the y values of all regular layers (0) and all buffer layers (1)
 * @param buffer_ranges a vector containing a tuple of the indices of the first and last node belonging to a certain buffer
 * @return sim_data_tuple see documentation of sim_data_tuple
 */
sim_data_tuple parallel_swap_framework::perform_swap_stream_and_collide_debug
(
    const std::vector<start_end_it_tuple> &fluid_nodes,
    const border_swap_information &bsi,
    std::vector<double> &distribution_values,    
    const access_function access_function,
    const std::tuple<std::vector<unsigned int>, std::vector<unsigned int>> &y_values,
    const std::vector<std::tuple<unsigned int, unsigned int>> &buffer_ranges
)
{
    std::cout << "Distribution values before stream and collide: " << std::endl;
    to_console::buffered::print_distribution_values(distribution_values, access_function);
    std::cout << std::endl;

    std::vector<velocity> velocities(TOTAL_NODE_COUNT, velocity{0,0});
    std::vector<double> densities(TOTAL_NODE_COUNT, -1);

    std::vector<double> corner_values = parallel_swap_framework::extract_corner_distributions(distribution_values, access_function);
    std::cout << "Extracted corner values " << std::endl;
    to_console::print_vector(corner_values, 5);

    // Border node initialization
    for(const auto node : bsi)
    {
        for(auto it = node.begin() + 1; it < node.end(); ++it)
        {
            swap_sequential::perform_swap_step(distribution_values, node[0], access_function, *it);
        }
    }
    std::cout << "Distribution values after border node initialization: " << std::endl;
    to_console::buffered::print_distribution_values(distribution_values, access_function);

    /* Buffer update */
    std::cout << "Copying to buffer" << std::endl;
    for(auto buffer_index = 0; buffer_index < BUFFER_COUNT; ++buffer_index)
    {
        parallel_swap_framework::swap_buffer_update(buffer_ranges[buffer_index], distribution_values, access_function);
    }
    std::cout << "Distribution values after buffer update: " << std::endl;
    to_console::buffered::print_distribution_values(distribution_values, access_function);

    // Swapping step
    for(auto subdomain = 0; subdomain < SUBDOMAIN_COUNT; ++subdomain)
    {
        for(auto node = std::get<0>(fluid_nodes[subdomain]); node <= std::get<1>(fluid_nodes[subdomain]); ++node)
        {
            swap_sequential::perform_swap_step(distribution_values, *node, access_function, swap_sequential::ACTIVE_STREAMING_DIRECTIONS);
        }
    }

    std::cout << "Distribution values after swap for every node: " << std::endl;
    to_console::buffered::print_distribution_values(distribution_values, access_function);
    std::cout << std::endl;

    // Restore precious order in here
    for(auto subdomain = 0; subdomain < SUBDOMAIN_COUNT; ++subdomain)
    {
        for(auto node = std::get<0>(fluid_nodes[subdomain]); node <= std::get<1>(fluid_nodes[subdomain]); ++node)
        {
            swap_sequential::restore_order(distribution_values, *node, access_function);
        }
    }

    std::cout << "Distribution values after ORDER has been restored for every node: " << std::endl;
    to_console::buffered::print_distribution_values(distribution_values, access_function);

    /* Perform collision for all fluid nodes */
    for(auto subdomain = 0; subdomain < SUBDOMAIN_COUNT; ++subdomain)
    {
        parallel_framework::perform_collision(fluid_nodes[subdomain], distribution_values, access_function, velocities, densities);
    }
    std::cout << "Distribution values after collision: " << std::endl;
    to_console::buffered::print_distribution_values(distribution_values, access_function);
    std::cout << std::endl;

    // Restore inout correctness
    parallel_swap_framework::restore_corner_distributions(corner_values, distribution_values, access_function);

    /* Update ghost nodes */
    parallel_framework::update_velocity_input_density_output(y_values, distribution_values, velocities, densities, access_function);
    std::cout << "Distribution values after ghost node update: " << std::endl;
    to_console::buffered::print_distribution_values(distribution_values, access_function);
    std::cout << std::endl;

    sim_data_tuple result{velocities, densities};

    return result;
}

/**
 * @brief Performs the streaming and collision step for all fluid nodes within the simulation domain.
 *        The border conditions are enforced through ghost nodes.
 *        This variant will print several debug comments to the console.
 * 
 * @param fluid_nodes a vector containing the first and last element of an iterator over all fluid nodes within each subdomain
 * @param bsi see documentation of border_swap_information
 * @param distribution_values a vector containing all distribution values
 * @param access_function the access to node values will be performed according to this access function
 * @param y_values a tuple containing the y values of all regular layers (0) and all buffer layers (1)
 * @param buffer_ranges a vector containing a tuple of the indices of the first and last node belonging to a certain buffer
 * @return sim_data_tuple see documentation of sim_data_tuple
 */
sim_data_tuple parallel_swap_framework::parallel_swap_stream_and_collide
(
    const std::vector<start_end_it_tuple> &fluid_nodes,
    const border_swap_information &bsi,
    std::vector<double> &distribution_values,    
    const access_function access_function,
    const std::tuple<std::vector<unsigned int>, std::vector<unsigned int>> &y_values,
    const std::vector<std::tuple<unsigned int, unsigned int>> &buffer_ranges
)
{
    std::vector<velocity> velocities(TOTAL_NODE_COUNT, velocity{0,0});
    std::vector<double> densities(TOTAL_NODE_COUNT, -1);
    std::vector<double> corner_values = parallel_swap_framework::extract_corner_distributions(distribution_values, access_function);

    /* Border node initialization */
    hpx::for_each
    (
        hpx::execution::par, 
        bsi.begin(), 
        bsi.end(), 
        [&](std::vector<unsigned int> node)
        {
        for(auto it = node.begin() + 1; it < node.end(); ++it)
            {
                swap_sequential::perform_swap_step(distribution_values, node[0], access_function, *it);
            }
        });

    /* Buffer update */
    hpx::experimental::for_loop(
        hpx::execution::par, 0, BUFFER_COUNT, 
        [&](unsigned int buffer_index)
        {
            parallel_swap_framework::swap_buffer_update(buffer_ranges[buffer_index], distribution_values, access_function);
        });

    hpx::experimental::for_loop(
        hpx::execution::par, 0, SUBDOMAIN_COUNT, 
        [&](unsigned int subdomain)
        {
            for(auto node = std::get<0>(fluid_nodes[subdomain]); node <= std::get<1>(fluid_nodes[subdomain]); ++node)
            {
                // Swapping step
                swap_sequential::perform_swap_step(distribution_values, *node, access_function, swap_sequential::ACTIVE_STREAMING_DIRECTIONS);

                // Restore precious order in here
                swap_sequential::restore_order(distribution_values, *node, access_function);

                /* Perform collision for all fluid nodes */
                parallel_framework::perform_collision(*node, distribution_values, access_function, velocities, densities);
            }
        });

    // for(auto subdomain = 0; subdomain < SUBDOMAIN_COUNT; ++subdomain)
    // {
    //     for(auto node = std::get<0>(fluid_nodes[subdomain]); node <= std::get<1>(fluid_nodes[subdomain]); ++node)
    //     {
    //         // Swapping step
    //         swap_sequential::perform_swap_step(distribution_values, *node, access_function, swap_sequential::ACTIVE_STREAMING_DIRECTIONS);

    //         // Restore precious order in here
    //         swap_sequential::restore_order(distribution_values, *node, access_function);

    //         /* Perform collision for all fluid nodes */
    //         parallel_framework::perform_collision(*node, distribution_values, access_function, velocities, densities);
    //     }
    // }

    // Restore inout correctness
    parallel_swap_framework::restore_corner_distributions(corner_values, distribution_values, access_function);

    /* Update ghost nodes */
    parallel_framework::update_velocity_input_density_output(y_values, distribution_values, velocities, densities, access_function);

    sim_data_tuple result{velocities, densities};

    return result;
}

void parallel_swap_framework::swap_buffer_update
(
    const std::tuple<unsigned int, unsigned int> &buffer_bounds,
    std::vector<double> &distribution_values,
    const access_function access_function
)
{
    unsigned int start = std::get<0>(buffer_bounds);
    unsigned int end = std::get<1>(buffer_bounds);

    // Clone values facin southward from subdomain above
    for(auto buffer_node = start; buffer_node <= end; ++buffer_node)
    {
        for(auto direction : {0,1,2})
        {
            distribution_values[access_function(buffer_node, direction)] = distribution_values[access_function(lbm_access::get_neighbor(buffer_node, 7), direction)];
        }
    }

    // Perform streaming across buffer
    for(auto buffer_node = start + 1; buffer_node <= end - 1; ++buffer_node)
    {
        for(auto direction : {6,7,8})
        {
            distribution_values[access_function(lbm_access::get_neighbor(buffer_node, direction), invert_direction(direction))] = distribution_values[access_function(lbm_access::get_neighbor(buffer_node, 1), direction)];
        }
    }
}

void parallel_swap_framework::restore_inout_correctness
(
    std::vector<double> &distribution_values,    
    const access_function access_function
)
{
    unsigned int restore_node = lbm_access::get_node_index(1,1);
    distribution_values[access_function(lbm_access::get_neighbor(restore_node, 0), 8)] = distribution_values[access_function(restore_node, 0)]; 
    restore_node = lbm_access::get_node_index(HORIZONTAL_NODES - 2, 1);
    distribution_values[access_function(lbm_access::get_neighbor(restore_node, 2), 6)] = distribution_values[access_function(restore_node, 2)]; 
    restore_node = lbm_access::get_node_index(1, VERTICAL_NODES - 2);
    distribution_values[access_function(lbm_access::get_neighbor(restore_node, 6), 2)] = distribution_values[access_function(restore_node, 6)]; 
    restore_node = lbm_access::get_node_index(HORIZONTAL_NODES - 2, VERTICAL_NODES - 2);
    distribution_values[access_function(lbm_access::get_neighbor(restore_node, 8), 0)] = distribution_values[access_function(restore_node, 8)]; 
}

std::vector<double> parallel_swap_framework::extract_corner_distributions
(
    const std::vector<double> &distribution_values,    
    const access_function access_function
)
{
    std::vector<double> result(4, 0);
    result[0] = distribution_values[access_function(lbm_access::get_node_index(0,0), 8)];
    result[1] = distribution_values[access_function(lbm_access::get_node_index(HORIZONTAL_NODES - 1, 0), 6)];
    result[2] = distribution_values[access_function(lbm_access::get_node_index(0, VERTICAL_NODES - 1), 2)];
    result[3] = distribution_values[access_function(lbm_access::get_node_index(HORIZONTAL_NODES - 1, VERTICAL_NODES - 1), 0)];
    return result;
}

void parallel_swap_framework::restore_corner_distributions
(
    const std::vector<double> &corner_values,  
    std::vector<double> &distribution_values,    
    const access_function access_function
)
{
    distribution_values[access_function(lbm_access::get_node_index(0,0), 8)] = corner_values[0];
    distribution_values[access_function(lbm_access::get_node_index(HORIZONTAL_NODES - 1, 0), 6)] = corner_values[1];
    distribution_values[access_function(lbm_access::get_node_index(0, VERTICAL_NODES - 1), 2)] = corner_values[2];
    distribution_values[access_function(lbm_access::get_node_index(HORIZONTAL_NODES - 1, VERTICAL_NODES - 1), 0)] = corner_values[3];
}