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
#include "../include/two_lattice_sequential.hpp"

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
void parallel_two_step_framework::run
(  
    const std::vector<start_end_it_tuple> &fluid_nodes,       
    std::vector<double> &distribution_values, 
    const border_swap_information &bsi,
    const access_function access_function,
    const unsigned int iterations
)
{
    to_console::print_run_greeting("parallel two-step algorithm", iterations);

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
        result[time] = parallel_two_step_framework::perform_ts_stream_and_collide_debug
        (fluid_nodes, bsi, distribution_values, access_function, y_values, buffer_ranges);
        std::cout << "\tFinished iteration " << time << std::endl;
    }
    to_console::buffered::print_simulation_results(result);
    std::cout << "All done, exiting simulation. " << std::endl;
}

/**
 * @brief Performs the streaming step for all fluid nodes within the specified bounds.
 * 
 * @param fluid_nodes a tuple of the first and last element of an iterator over all fluid nodes within the respective subdomain
 * @param distribution_values a vector containing all distribution distribution_values
 * @param access_function the access to node values will be performed according to this access function
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
 * @brief Performs the collision step for all fluid nodes within the specified bounds.
 * 
 * @param fluid_node_bounds a tuple of the first and last element of an iterator over all fluid nodes within the respective subdomain
 * @param distribution_values a vector containing all distribution distribution_values
 * @param access_function the access to node values will be performed according to this access function
 * @param velocities a vector containing the velocity values of all nodes
 * @param densities a vector containing the density values of all nodes
 */
void parallel_two_step_framework::perform_collision
(
    const start_end_it_tuple fluid_node_bounds,
    std::vector<double> &distribution_values, 
    const access_function &access_function, 
    std::vector<velocity> &velocities, 
    std::vector<double> &densities
)
{
    std::vector<double> current_distributions(DIRECTION_COUNT, 0);
    velocity current_velocity = {0,0};
    double current_density = 0;
    for(auto it = std::get<0>(fluid_node_bounds); it <= std::get<1>(fluid_node_bounds); ++it)
    {
        if(*it == 82)
        {
            std::cout << std::setprecision(8) << std::fixed;
            current_distributions = 
                lbm_access::get_distribution_values_of(distribution_values, *it, access_function);
            std::cout << "Two-step: Performing collision for node 82 " << std::endl;
            std::cout << "Got distribution values " << std::endl;
            to_console::print_vector(lbm_access::get_distribution_values_of(distribution_values, *it, access_function), DIRECTION_COUNT);
            std::cout << "Got velocity (" << std::get<0>(macroscopic::flow_velocity(current_distributions)) << ", " << std::get<1>(macroscopic::flow_velocity(current_distributions)) << ")" << std::endl;
            std::cout << "Got density: " << macroscopic::density(current_distributions) << std::endl;
            std::cout << "Resulting distribution is " << std::endl;
            to_console::print_vector(collision::collide_bgk(current_distributions, current_velocity, current_density), DIRECTION_COUNT);
            std::cout << std::setprecision(3) << std::fixed;
        }
        current_distributions = 
            lbm_access::get_distribution_values_of(distribution_values, *it, access_function);
        current_velocity = macroscopic::flow_velocity(current_distributions);    
        velocities[*it] = current_velocity;
        current_density = macroscopic::density(current_distributions);
        densities[*it] = current_density;
        current_distributions = collision::collide_bgk(current_distributions, current_velocity, current_density);
        lbm_access::set_distribution_values_of(current_distributions, distribution_values, *it, access_function);
    }
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
sim_data_tuple parallel_two_step_framework::perform_ts_stream_and_collide_debug
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

    /* Streaming */
    for(auto subdomain = 0; subdomain < SUBDOMAIN_COUNT; ++subdomain)
    {
        parallel_two_step_framework::perform_stream(fluid_nodes[subdomain], distribution_values, access_function);
    }

    /* Get remaining streams from buffer */
    std::cout << "Copying from buffer" << std::endl;

    for(auto buffer_index = 0; buffer_index < BUFFER_COUNT; ++buffer_index)
    {
        parallel_framework::copy_from_buffer(std::make_tuple(std::get<0>(buffer_ranges[buffer_index])+1, std::get<1>(buffer_ranges[buffer_index])-1), distribution_values, access_function);
    }
    std::cout << "\t Distribution values after streaming:" << std::endl;
    to_console::buffered::print_distribution_values(distribution_values, access_function);
    std::cout << std::endl;

    /* Perform bounce-back using ghost nodes */
    parallel_two_step_framework::perform_boundary_update(bsi, distribution_values, access_function);
    std::cout << "\t Distribution values after bounce-back update:" << std::endl;
    to_console::buffered::print_distribution_values(distribution_values, access_function);
    std::cout << std::endl;

    /* Perform inflow and outflow using ghost nodes */
    std::cout << "Performing ghost stream inout" << std::endl;
    parallel_two_step_framework::ghost_stream_inout(distribution_values, access_function, y_values);
    
    std::cout << "\t Distribution values after inflow and outflow via ghost nodes:" << std::endl;
    to_console::buffered::print_distribution_values(distribution_values, access_function);
    std::cout << std::endl;

    double current_density = 0;
    velocity current_velocity = {0,0};
    std::vector<double> current_distributions(DIRECTION_COUNT, 0);

    /* Perform collision for all fluid nodes */
    for(auto subdomain = 0; subdomain < SUBDOMAIN_COUNT; ++subdomain)
    {
        //parallel_two_step_framework::perform_collision(fluid_nodes[subdomain], distribution_values, access_function, velocities, densities);
        for(auto it = std::get<0>(fluid_nodes[subdomain]); it <= std::get<1>(fluid_nodes[subdomain]); ++it)
        {   
            current_distributions = 
                lbm_access::get_distribution_values_of(distribution_values, *it, access_function);

            current_velocity = macroscopic::flow_velocity(current_distributions);    
            velocities[*it] = current_velocity;

            current_density = macroscopic::density(current_distributions);
            densities[*it] = current_density;
            
            two_lattice_sequential::tl_collision(
                distribution_values, 
                *it, 
                current_distributions,
                access_function, 
                current_velocity, 
                current_density);
        }
    }
    std::cout << "Distribution values after collision: " << std::endl;
    to_console::buffered::print_distribution_values(distribution_values, access_function);
    std::cout << std::endl;

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
 * 
 * @param fluid_nodes a vector containing the first and last element of an iterator over all fluid nodes within each subdomain
 * @param bsi see documentation of border_swap_information
 * @param distribution_values a vector containing all distribution values
 * @param access_function the access to node values will be performed according to this access function
 * @param y_values a tuple containing the y values of all regular layers (0) and all buffer layers (1)
 * @param buffer_ranges a vector containing a tuple of the indices of the first and last node belonging to a certain buffer
 * @return sim_data_tuple see documentation of sim_data_tuple
 */
sim_data_tuple parallel_two_step_framework::parallel_ts_stream_and_collide
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

    /* Perform streaming for all fluid nodes */
    hpx::experimental::for_loop(
        hpx::execution::par, 0, SUBDOMAIN_COUNT,
        [&](int subdomain)
        {  
            parallel_two_step_framework::perform_stream(fluid_nodes[subdomain], distribution_values, access_function);
        });

    /* Get remaining streams from buffer */
    hpx::experimental::for_loop(
        hpx::execution::par, 0, BUFFER_COUNT, 
        [&](unsigned int buffer_index)
        {
            parallel_framework::copy_from_buffer(buffer_ranges[buffer_index], distribution_values, access_function);
        }
    );

    /* Perform bounce-back using ghost nodes */
    parallel_two_step_framework::perform_boundary_update(bsi, distribution_values, access_function);

    /* Perform inflow and outflow using ghost nodes */
    parallel_two_step_framework::ghost_stream_inout(distribution_values, access_function, y_values);

    /* Perform collision for all fluid nodes */
    hpx::experimental::for_loop(
        hpx::execution::par, 0, SUBDOMAIN_COUNT,
        [&](int subdomain)
        {
        
        parallel_framework::perform_collision(fluid_nodes[subdomain], distribution_values, access_function, velocities, densities);
    });

    /* Update ghost nodes */
    parallel_framework::update_velocity_input_density_output(y_values, distribution_values, velocities, densities, access_function);

    sim_data_tuple result{velocities, densities};

    return result;
}

/**
 * @brief Realizes inflow and outflow by an inward stream of each border node.
 * 
 * @param distribution_values a vector containing the distribution values of all nodes
 * @param access_function the access function used to access the distribution values
 * @param y_values a tuple containing the y values of all regular layers (0) and all buffer layers (1)
 */
void parallel_two_step_framework::ghost_stream_inout
(
    std::vector<double> &distribution_values, 
    const access_function access_function,
    const std::tuple<std::vector<unsigned int>, std::vector<unsigned int>> &y_values
)
{

    // // Restore buffer correctness
    // hpx::for_each
    // (
    //     hpx::execution::par, 
    //     std::get<1>(y_values).begin(), 
    //     std::get<1>(y_values).end(), 
    //     [&](int y)
    //     {
    //         unsigned int current_border_node = lbm_access::get_node_index(0,y);
    //         parallel_framework::copy_to_buffer_node(current_border_node, distribution_values, access_function);
        
    //         current_border_node = lbm_access::get_node_index(HORIZONTAL_NODES - 1,y);
    //         parallel_framework::copy_to_buffer_node(current_border_node, distribution_values, access_function);

    //         // unsigned int buffer_node = 0; 
    //         // unsigned int current_neighbor = 0;

    //         // for(auto side : {1, HORIZONTAL_NODES - 2})
    //         // {
    //         //     buffer_node = lbm_access::get_node_index(side, y);
    //         //     current_neighbor = lbm_access::get_neighbor(buffer_node, 1);
    //         //     distribution_values[access_function(current_neighbor, 2)] = distribution_values[access_function(buffer_node, 2)];

    //         //     current_neighbor = lbm_access::get_neighbor(buffer_node, 7);
    //         //     distribution_values[access_function(current_neighbor, 8)] = distribution_values[access_function(buffer_node, 8)];
    //         // }
    //     }
    // );

    hpx::for_each
    (
        hpx::execution::par, 
        std::get<0>(y_values).begin(), 
        std::get<0>(y_values).end(), 
        [&](int y)
        {
            // Update inlets
            unsigned int current_border_node = lbm_access::get_node_index(1,y);
            for(const auto direction : inflow_instream_dirs)
            {
                distribution_values[access_function(current_border_node, direction)] = 
                    distribution_values[access_function(lbm_access::get_neighbor(current_border_node, invert_direction(direction)), direction)];
            }

            // Update outlets
            current_border_node = lbm_access::get_node_index(HORIZONTAL_NODES - 2,y);
            for(const auto direction : outflow_instream_dirs)
            {
                distribution_values[access_function(current_border_node, direction)] = 
                    distribution_values[access_function(lbm_access::get_neighbor(current_border_node, invert_direction(direction)), direction)];
            }
        }
    );


}

/**
 * @brief Performs a halfway bounce-back streaming update for all fluid nodes within the simulation domain.
 * 
 * @param bsi see documentation of border_swap_information
 * @param distribution_values a vector containing the distribution values of all nodes
 * @param access_function the access function used to access the distribution values
 */
void parallel_two_step_framework::perform_boundary_update
(
    const border_swap_information &bsi,
    std::vector<double> &distribution_values, 
    const access_function access_function
)
{
    // hpx::for_each
    // (
    //     hpx::execution::par, 
    //     bsi.begin(), 
    //     bsi.end(), 
    //     [&](std::vector<unsigned int> current)
    //     {
    //         int current_border_node = current[0];
    //         std::set<unsigned int> remaining_dirs = bounce_back::determine_bounce_back_directions(current);
    //         std::cout << "Performing bounce-back for node " << current_border_node << " in directions ";
    //         to_console::print_set(remaining_dirs);
    //         for(const auto direction : remaining_dirs)
    //         {
    //             distribution_values[access_function(current_border_node, direction)] = 
    //             distribution_values[access_function(lbm_access::get_neighbor(current_border_node, invert_direction(direction)), invert_direction(direction))];
    //         }
    //     }
    // );
    
    for(auto current : bsi)
    {
        int current_border_node = current[0];
        std::set<unsigned int> remaining_dirs = bounce_back::determine_bounce_back_directions(current);
        std::cout << "Performing bounce-back for node " << current_border_node << " in directions ";
        to_console::print_set(remaining_dirs);
        for(const auto direction : remaining_dirs)
        {
            distribution_values[access_function(current_border_node, direction)] = 
            distribution_values[access_function(lbm_access::get_neighbor(current_border_node, invert_direction(direction)), invert_direction(direction))];
        }
    }



}
