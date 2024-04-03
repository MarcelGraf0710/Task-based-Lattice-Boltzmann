#include "../include/parallel_two_lattice_framework.hpp"

#include "../include/two_lattice_sequential.hpp"
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


    // hpx::for_each
    // (
    //     hpx::execution::par, 
    //     fluid_nodes.begin(), 
    //     fluid_nodes.end(), 
    //     [&source, &destination, access_function, &velocities, &densities](unsigned int fluid_node)
    //     {
    //         tl_stream_and_collide_helper(source, destination, access_function, fluid_node, velocities, densities);
    //     }
    // );

/**
 * @brief Performs the sequential two-lattice algorithm for the specified number of iterations.
 * 
 * @param fluid_nodes A vector containing the indices of all fluid nodes in the domain
 * @param boundary_nodes see documentation of border_swap_information
 * @param distribution_values_0 source for even time steps and destination for odd time steps
 * @param distribution_values_1 source for odd time steps and destination for even time steps
 * @param access_function the access function according to which distribution values are to be accessed
 * @param iterations this many iterations will be performed
 */
void parallel_two_lattice_framework::run
(  
    std::vector<start_end_it_tuple> &fluid_nodes,       
    const border_swap_information &boundary_nodes,
    std::vector<double> &distribution_values_0, 
    std::vector<double> &distribution_values_1,   
    const access_function access_function,
    const unsigned int iterations
)
{
    to_console::print_run_greeting("Parallel two-lattice algorithm (framework version)", iterations);

    std::vector<double> &source = distribution_values_0;
    std::vector<double> &destination = distribution_values_1;
    std::vector<double> &temp = distribution_values_1;

    std::vector<sim_data_tuple>result(
        iterations, 
        std::make_tuple(std::vector<velocity>(TOTAL_NODE_COUNT, {0,0}), std::vector<double>(TOTAL_NODE_COUNT, 0)));

    for(auto time = 0; time < iterations; ++time)
    {
        std::cout << "\033[33mIteration " << time << ":\033[0m" << std::endl;
        // Global boundary update
        bounce_back::emplace_bounce_back_values(boundary_nodes, source, access_function);

        // Buffer update
        for (auto buffer_index = 0; buffer_index < BUFFER_COUNT; ++buffer_index)
        {
            parallel_framework::buffer_copy_update(parallel_framework::get_buffer_node_range(buffer_index), source, access_function);
        }

        // TODO: Implement stream and collide
        result[time] = parallel_two_lattice_framework::perform_tl_stream_and_collide_parallel
        (
            fluid_nodes, boundary_nodes, source, destination, access_function
        );

        std::cout << "\tFinished iteration " << time << std::endl;
        std::swap(source, destination);
    }

    to_console::print_simulation_results_buffered(result);
    std::cout << "All done, exiting simulation. " << std::endl;
}


/**
 * @brief Performs the combined streaming and collision step for all fluid nodes within the simulation domain.
 *        The border conditions are enforced through ghost nodes.
 *        This variant of the combined streaming and collision step will print several debug comments to the console.
 * 
 * @param fluid_nodes a vector containing the indices of all fluid nodes within the simulation domain.
 * @param bsi see documentation of border_swap_information
 * @param source a vector containing the distribution values of the previous time step
 * @param destination the distribution values will be written to this vector after performing both steps.
 * @param access_function the function used to access the distribution values
 * @return see documentation of sim_data_tuple
 */
sim_data_tuple parallel_two_lattice_framework::perform_tl_stream_and_collide_debug
(
    std::vector<start_end_it_tuple> &fluid_nodes,
    const border_swap_information &bsi,
    std::vector<double> &source, 
    std::vector<double> &destination,    
    const access_function access_function
)
{
    std::vector<velocity> velocities(TOTAL_NODE_COUNT, velocity{0,0});
    velocity current_velocity = {0,0};

    std::vector<double> densities(TOTAL_NODE_COUNT, -1);
    double current_density = 0;

    std::vector<double> current_distributions(DIRECTION_COUNT, 0);

    start_end_it_tuple bounds;

    std::cout << "\t TL stream and collide: initializations and declarations performed." << std::endl;

    for (auto subdomain = 0; subdomain < SUBDOMAIN_COUNT; ++subdomain)
    {
        std::cout << "\t\033[33mPerforming iteration for subdomain " << subdomain << "\033[0m" << std::endl;
        std::cout << std::endl;

        bounds = fluid_nodes[subdomain];

        std::cout << "\t SOURCE before stream and collide: " << std::endl;
        to_console::print_distribution_values_buffered(source, access_function);

        std::cout << "DESTINATION as received by perform_tl_stream_and_collide: " << std::endl;
        to_console::print_distribution_values_buffered(destination, access_function);

        /* Streaming step */
        for(auto it = std::get<0>(bounds); it <= std::get<1>(bounds); ++it)
        {
            two_lattice_sequential::tl_stream(
                source, 
                destination, 
                access_function, 
                *it);
        }
        std::cout << "DESTINATION after streaming: " << std::endl;
        to_console::print_distribution_values_buffered(destination, access_function);

        /* Collision step */
        for(auto it = std::get<0>(bounds); it <= std::get<1>(bounds); ++it)
        {   
            current_distributions = 
                lbm_access::get_distribution_values_of(destination, *it, access_function);

            current_velocity = macroscopic::flow_velocity(current_distributions);    
            velocities[*it] = current_velocity;

            current_density = macroscopic::density(current_distributions);
            densities[*it] = current_density;
            
            two_lattice_sequential::tl_collision(
                destination, 
                *it, 
                current_distributions,
                access_function, 
                current_velocity, 
                current_density);
        }
        std::cout << "\t DESTINATION after collision: " << std::endl;
        to_console::print_distribution_values_buffered(destination, access_function);
    }

    boundary_conditions::update_velocity_input_density_output(destination, velocities, densities, access_function);
    std::cout << "Updated inlet and outlet ghost nodes." <<std::endl;
    to_console::print_distribution_values_buffered(destination, access_function);

    sim_data_tuple result{velocities, densities};

    return result;
}

/**
 * @brief Performs the combined streaming and collision step for all fluid nodes within the simulation domain.
 *        The border conditions are enforced through ghost nodes.
 *        This variant of the combined streaming and collision step will print several debug comments to the console.
 * 
 * @param fluid_nodes a vector containing the indices of all fluid nodes within the simulation domain.
 * @param bsi see documentation of border_swap_information
 * @param source a vector containing the distribution values of the previous time step
 * @param destination the distribution values will be written to this vector after performing both steps.
 * @param access_function the function used to access the distribution values
 * @return see documentation of sim_data_tuple
 */
sim_data_tuple parallel_two_lattice_framework::perform_tl_stream_and_collide_parallel
(
    std::vector<start_end_it_tuple> &fluid_nodes,
    const border_swap_information &bsi,
    std::vector<double> &source, 
    std::vector<double> &destination,    
    const access_function access_function
)
{
    std::vector<velocity> velocities(TOTAL_NODE_COUNT, velocity{0,0});


    std::vector<double> densities(TOTAL_NODE_COUNT, -1);

    

    start_end_it_tuple bounds;

    std::vector<int> subdomains(SUBDOMAIN_COUNT, 0);
    for(auto i = 0; i < SUBDOMAIN_COUNT; ++i)
    {
        subdomains[i] = i;
    }

    hpx::for_each
    (
        hpx::execution::par, 
        subdomains.begin(), 
        subdomains.end(), 
        [&source, &destination, access_function, &velocities, &densities, &bounds, &fluid_nodes](unsigned int fluid_node)
        {
            bounds = fluid_nodes[fluid_node];
            parallel_two_lattice_framework::tl_stream_and_collide_helper(source, destination, access_function, fluid_node, velocities, densities, bounds);
        }
    );

    boundary_conditions::update_velocity_input_density_output(destination, velocities, densities, access_function);

    sim_data_tuple result{velocities, densities};

    return result;
}


void parallel_two_lattice_framework::tl_stream_and_collide_helper
(
    std::vector<double> &source, 
    std::vector<double> &destination, 
    const access_function &access_function, 
    const unsigned int fluid_node, 
    std::vector<velocity> &velocities, 
    std::vector<double> &densities,
    start_end_it_tuple bounds
)
{
    std::vector<double> current_distributions(DIRECTION_COUNT, 0);
    velocity current_velocity = {0,0};
    double current_density = 0;

    for(auto it = std::get<0>(bounds); it <= std::get<1>(bounds); ++it)
    {
        /* Streaming step */
        two_lattice_sequential::tl_stream(
            source, 
            destination, 
            access_function, 
            *it);

        /* Collision step */
        current_distributions = 
            lbm_access::get_distribution_values_of(destination, *it, access_function);

        current_velocity = macroscopic::flow_velocity(current_distributions);    
        velocities[*it] = current_velocity;

        current_density = macroscopic::density(current_distributions);
        densities[*it] = current_density;
        
        two_lattice_sequential::tl_collision(
            destination, 
            *it, 
            current_distributions,
            access_function, 
            current_velocity, 
            current_density);           
    }
} 