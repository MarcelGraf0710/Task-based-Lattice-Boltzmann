#include "../include/two_step_sequential.hpp"
#include "../include/access.hpp"
#include "../include/boundaries.hpp"
#include "../include/macroscopic.hpp"
#include "../include/collision.hpp"
#include <iostream>

/**
 * @brief Performs the streaming step for all fluid nodes within the simulation domain.
 *        Notice that each node is streaming outward.
 * 
 * @param fluid_nodes A vector containing the indices of all fluid nodes in the domain
 * @param values a vector containing all distribution values
 * @param access_function the access to node values will be performed according to this access function.
 */
void two_step_sequential::perform_fast_stream
(
    const std::vector<unsigned int> &fluid_nodes, 
    std::vector<double> &values, 
    const access_function access_function
)
{
    // All directions that require left-to-right and/or bottom-to-top node iteration order
    for(auto it = fluid_nodes.begin(); it < fluid_nodes.end(); ++it)
    {
        values[access_function(access::get_neighbor(*it, 0), 0)] = values[access_function(*it, 0)];
        values[access_function(access::get_neighbor(*it, 1), 1)] = values[access_function(*it, 1)];
        values[access_function(access::get_neighbor(*it, 2), 2)] = values[access_function(*it, 2)];
        values[access_function(access::get_neighbor(*it, 3), 3)] = values[access_function(*it, 3)];
    }

    // All directions that require right-to-left and/or top-to-bottom node iteration order
    for(auto it = fluid_nodes.end() - 1; it >= fluid_nodes.begin(); --it)
    {
        values[access_function(access::get_neighbor(*it, 5), 5)] = values[access_function(*it, 5)];
        values[access_function(access::get_neighbor(*it, 6), 6)] = values[access_function(*it, 6)];
        values[access_function(access::get_neighbor(*it, 7), 7)] = values[access_function(*it, 7)];
        values[access_function(access::get_neighbor(*it, 8), 8)] = values[access_function(*it, 8)];
    }
}

/**
 * @brief Performs the combined streaming and collision step for all fluid nodes within the simulation domain.
 *        The border conditions are enforced through ghost nodes.
 *        This variant of the combined streaming and collision step will print several debug comments to the console.
 */
sim_data_tuple two_step_sequential::perform_ts_stream_and_collide_debug
(
    const std::vector<unsigned int> &fluid_nodes,
    const border_swap_information &bsi,
    std::vector<double> &distribution_values,    
    const access_function access_function
)
{
    std::cout << "Distribution values before stream and collide: " << std::endl;
    to_console::print_distribution_values(distribution_values, access_function);
    std::cout << std::endl;

    std::vector<velocity> velocities(TOTAL_NODE_COUNT, velocity{0,0});
    std::vector<double> densities(TOTAL_NODE_COUNT, -1);
    std::vector<double> current_distributions(DIRECTION_COUNT, 0);

    /* Streaming */
    two_step_sequential::perform_fast_stream(fluid_nodes, distribution_values, access_function);
    std::cout << "\t Distribution values after streaming:" << std::endl;
    to_console::print_distribution_values(distribution_values, access_function);
    std::cout << std::endl;

    /* Perform bounce-back using ghost nodes */
    bounce_back::perform_boundary_update(bsi, distribution_values, access_function);
    std::cout << "\t Distribution values after bounce-back update:" << std::endl;
    to_console::print_distribution_values(distribution_values, access_function);
    std::cout << std::endl;

    /* Perform inflow and outflow using ghost nodes */
    std::cout << "Performing ghost stream inout " << std::endl;
    boundary_conditions::ghost_stream_inout(distribution_values, access_function);
    std::cout << "\t Distribution values after inflow and outflow via ghost nodes:" << std::endl;
    to_console::print_distribution_values(distribution_values, access_function);
    std::cout << std::endl;

    /* Perform collision for all fluid nodes */
    velocities = macroscopic::calculate_all_velocities(fluid_nodes, distribution_values, access_function);
    densities = macroscopic::calculate_all_densities(fluid_nodes, distribution_values, access_function);
    collision::collide_all_bgk(fluid_nodes, distribution_values, velocities, densities, access_function);
    std::cout << "\t Distribution values after collision:" << std::endl;
    to_console::print_distribution_values(distribution_values, access_function);
    std::cout << std::endl;

    /* Update ghost nodes */
    boundary_conditions::update_density_input_density_output(distribution_values, access_function);
    std::cout << "Distribution values after ghost node update: " << std::endl;
    to_console::print_distribution_values(distribution_values, access_function);
    std::cout << std::endl;

    unsigned int update_node = 0;
    for(auto x = 0; x < HORIZONTAL_NODES; x = x + HORIZONTAL_NODES - 1)
    {
        for(auto y = 1; y < VERTICAL_NODES - 1; ++y)
        {
            update_node = access::get_node_index(x,y);
            current_distributions = 
                access::get_distribution_values_of(distribution_values, update_node, access_function);
                
            velocities[update_node] = macroscopic::flow_velocity(current_distributions);
            densities[update_node] = macroscopic::density(current_distributions);
        }
    }

    sim_data_tuple result{velocities, densities};

    return result;
}

/**
 * @brief Performs the sequential two-step algorithm for the specified number of iterations.
 * 
 * @param fluid_nodes A vector containing the indices of all fluid nodes in the domain
 * @param values the vector containing the distribution values of all nodes
 * @param access_function the access function according to which the values are to be accessed
 * @param iterations this many iterations will be performed
 */
void two_step_sequential::run
(  
    std::vector<unsigned int> &fluid_nodes,       
    std::vector<double> &values, 
    border_swap_information &bsi,
    access_function access_function,
    unsigned int iterations
)
{
    std::cout << "------------------------------------------------------------------------------------------------------------------------" << std::endl;
    std::cout << "Now running sequential two step algorithm for " << iterations << " iterations." << std::endl;
    std::cout << std::endl;

    std::vector<sim_data_tuple>result(
        iterations, 
        std::make_tuple(std::vector<velocity>(TOTAL_NODE_COUNT, {0,0}), std::vector<double>(TOTAL_NODE_COUNT, 0)));

    for(auto time = 0; time < iterations; ++time)
    {
        std::cout << "\033[33mIteration " << time << ":\033[0m" << std::endl;
        result[time] = two_step_sequential::perform_ts_stream_and_collide_debug
        (
            fluid_nodes, 
            bsi, 
            values, 
            access_function
        );
        std::cout << "Finished iteration " << time << std::endl;
    }

    std::cout << std::endl;
    std::cout << std::endl;

    std::cout << "Velocity values: " << std::endl;
    std::cout << std::endl;
    for(auto i = 0; i < iterations; ++i)
    {
        std::cout << "t = " << i << std::endl;
        std::cout << "-------------------------------------------------------------------------------- " << std::endl;
        to_console::print_velocity_vector(std::get<0>(result[i]));
        std::cout << std::endl;
    }
    std::cout << std::endl;
    std::cout << std::endl;

    std::cout << "Density values: " << std::endl;
    std::cout << std::endl;
    
    for(auto i = 0; i < iterations; ++i)
    {
        std::cout << "t = " << i << std::endl;
        std::cout << "-------------------------------------------------------------------------------- " << std::endl;
        to_console::print_vector(std::get<1>(result[i]));
        std::cout << std::endl;
    }
    std::cout << "All done, exiting simulation. " << std::endl;
   
}
