#include "../include/two_lattice_parallel.hpp"
#include "../include/access.hpp"
#include "../include/defines.hpp"
#include "../include/boundaries.hpp"
#include "../include/collision.hpp"
#include "../include/macroscopic.hpp"
#include "../include/utils.hpp"
#include <hpx/format.hpp>
#include <hpx/future.hpp>
#include <hpx/algorithm.hpp>
#include <hpx/execution.hpp>
#include <hpx/iostream.hpp>
#include <set>
#include <iostream>
#include <functional>

/**
 * @brief Performs the combined streaming and collision step for all fluid nodes within the simulation domain.
 *        The border conditions are enforced through ghost nodes.
 * 
 * @param fluid_nodes a vector containing the indices of all fluid nodes within the simulation domain.
 * @param bsi see documentation of border_swap_information
 * @param source a vector containing the distribution values of the previous time step
 * @param destination the distribution values will be written to this vector after performing both steps.
 * @param access_function the function used to access the distribution values
 * @return see documentation of sim_data_tuple
 */
sim_data_tuple two_lattice_parallel::perform_tl_stream_and_collide
(
    const std::vector<unsigned int> &fluid_nodes,
    const border_swap_information &bsi,
    std::vector<double> &source, 
    std::vector<double> &destination,    
    const access_function access_function
)
{
    std::vector<velocity> velocities(TOTAL_NODE_COUNT, velocity{0,0});
    std::vector<double> densities(TOTAL_NODE_COUNT, -1);

    /* Boundary node treatment */
    bounce_back::emplace_bounce_back_values(bsi, source, access_function);

    /* Combined stream and collision step */
    hpx::for_each
    (
        hpx::execution::par, 
        fluid_nodes.begin(), 
        fluid_nodes.end(), 
        std::bind(tl_stream_and_collide_helper, source, destination, access_function, std::placeholders::_1, velocities, densities)
    );

    boundary_conditions::update_density_input_density_output(destination, velocities, densities, access_function);

    sim_data_tuple result{velocities, densities};

    return result;
}

void two_lattice_parallel::tl_stream_and_collide_helper
(
    std::vector<double> &source, 
    std::vector<double> &destination, 
    const access_function &access_function, 
    const unsigned int fluid_node, 
    std::vector<velocity> &velocities, 
    std::vector<double> &densities
)
{
    char const* msg = "Currently treating node {1} on locality {2}\n";

    hpx::util::format_to(hpx::cout, msg, fluid_node, hpx::get_locality_id())
            << std::flush;
    //hpx::cout << "Currently treating node " << fluid_node << "\n" <<  std::flush;
    velocity current_velocity = {0,0};
    double current_density = 0;
    std::vector<double> current_distributions(DIRECTION_COUNT, 0);

    two_lattice_parallel::tl_stream(
        source,
        destination,
        access_function,
        fluid_node);

    current_distributions =
        lbm_access::get_distribution_values_of(destination, fluid_node, access_function);

    current_velocity = macroscopic::flow_velocity(current_distributions);
    velocities[fluid_node] = current_velocity;

    current_density = macroscopic::density(current_distributions);
    densities[fluid_node] = current_density;

    two_lattice_parallel::tl_collision(
        destination,
        fluid_node,
        current_distributions,
        access_function,
        current_velocity,
        current_density);
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
sim_data_tuple two_lattice_parallel::perform_tl_stream_and_collide_debug
(
    const std::vector<unsigned int> &fluid_nodes,
    const border_swap_information &bsi,
    std::vector<double> &source, 
    std::vector<double> &destination,    
    const access_function access_function
)
{
    std::cout << "\t SOURCE before stream and collide: " << std::endl;
    to_console::print_distribution_values(source, access_function);

    std::cout << "DESTINATION as received by perform_tl_stream_and_collide: " << std::endl;
    to_console::print_distribution_values(destination, access_function);

    std::vector<velocity> velocities(TOTAL_NODE_COUNT, velocity{0,0});
    //velocity current_velocity = {0,0};

    std::vector<double> densities(TOTAL_NODE_COUNT, -1);
    //double current_density = 0;

    //std::vector<double> current_distributions(DIRECTION_COUNT, 0);

    std::cout << "\t TL stream and collide: initializations and declarations performed." << std::endl;

    /* Boundary node treatment */
    bounce_back::emplace_bounce_back_values(bsi, source, access_function);
    std::cout << "SOURCE after emplace bounce-back values: " << std::endl;
    to_console::print_distribution_values(source, access_function);

    /* Streaming step I AM WORKING */
    hpx::for_each
    (
        hpx::execution::par, 
        fluid_nodes.begin(), 
        fluid_nodes.end(), 
        [&source, &destination, access_function](unsigned int fluid_node)
        {
            two_lattice_parallel::tl_stream(
            source, 
            destination, 
            access_function, 
            fluid_node);
        }
    );

    // for(const auto fluid_node : fluid_nodes)
    // {
    //     two_lattice_parallel::tl_stream(
    //         source, 
    //         destination, 
    //         access_function, 
    //         fluid_node);
    // }

    std::cout << "DESTINATION after streaming: " << std::endl;
    to_console::print_distribution_values(destination, access_function);

    /* Collision step */
    hpx::for_each
    (
        hpx::execution::seq, 
        fluid_nodes.begin(), 
        fluid_nodes.end(), 
        [/*&current_distributions, &current_velocity, &current_density, */ &velocities, &densities, &source, &destination, access_function](unsigned int fluid_node)
        {
            std::vector<double> current_distributions(DIRECTION_COUNT, 0);
            current_distributions = 
                lbm_access::get_distribution_values_of(destination, fluid_node, access_function);

            std::cout << "Current distribution values: \n" << std::endl;
            to_console::print_vector(current_distributions);

            velocity current_velocity = macroscopic::flow_velocity(current_distributions);    
            // std::cout << "Current velocity: \n" << std::flush;
            // std::cout << "(" << current_velocity[0] << ", " << current_velocity[1] << ")" << std::flush;
            velocities[fluid_node] = current_velocity;

            double current_density = macroscopic::density(current_distributions);
            // std::cout << "Current density: " << current_density << std::flush;
            densities[fluid_node] = current_density;
            
            two_lattice_parallel::tl_collision(
                destination, 
                fluid_node, 
                current_distributions,
                access_function, 
                current_velocity, 
                current_density);
        }
    );


    // for(const auto fluid_node : fluid_nodes)
    // {   
    //     current_distributions = 
    //         lbm_access::get_distribution_values_of(destination, fluid_node, access_function);

    //     current_velocity = macroscopic::flow_velocity(current_distributions);    
    //     velocities[fluid_node] = current_velocity;

    //     current_density = macroscopic::density(current_distributions);
    //     densities[fluid_node] = current_density;
        
    //     two_lattice_parallel::tl_collision(
    //         destination, 
    //         fluid_node, 
    //         current_distributions,
    //         access_function, 
    //         current_velocity, 
    //         current_density);
    // }

    std::cout << "\t DESTINATION after collision: " << std::endl;
    to_console::print_distribution_values(destination, access_function);

    boundary_conditions::update_velocity_input_density_output(destination, velocities, densities, access_function);
    std::cout << "Updated inlet and outlet ghost nodes." <<std::endl;
    to_console::print_distribution_values(destination, access_function);

    sim_data_tuple result{velocities, densities};

    return result;
}

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
void two_lattice_parallel::run
(  
    const std::vector<unsigned int> &fluid_nodes,       
    const border_swap_information &boundary_nodes,
    std::vector<double> &distribution_values_0, 
    std::vector<double> &distribution_values_1,   
    const access_function access_function,
    const unsigned int iterations
)
{
    to_console::print_run_greeting("sequential two-lattice algorithm", iterations);

    std::vector<double> &source = distribution_values_0;
    std::vector<double> &destination = distribution_values_1;
    std::vector<double> &temp = distribution_values_1;
    std::vector<sim_data_tuple>result(
        iterations, 
        std::make_tuple(std::vector<velocity>(TOTAL_NODE_COUNT, {0,0}), std::vector<double>(TOTAL_NODE_COUNT, 0)));

    for(auto time = 0; time < iterations; ++time)
    {
        std::cout << "\033[33mIteration " << time << ":\033[0m" << std::endl;
        result[time] = two_lattice_parallel::perform_tl_stream_and_collide_debug
        (
            fluid_nodes, 
            boundary_nodes, 
            source, 
            destination, 
            access_function
        );     
        std::cout << "\tFinished iteration " << time << std::endl;
        std::swap(source, destination);
    }

    to_console::print_simulation_results(result);
    std::cout << "All done, exiting simulation. " << std::endl;
}

/**
 * @brief Determines the remaining streaming option for a node based on the specified border 
 *        information vector.
 * 
 * @param current_border_info an entry of a border_swap_information object
 * @return a set containing all remaining streaming directions
 */
std::set<unsigned int> two_lattice_parallel::determine_streaming_directions
(
    const std::vector<unsigned int> &current_border_info
)
{
    std::set<unsigned int> remaining_dirs = {STREAMING_DIRECTIONS.begin(), STREAMING_DIRECTIONS.end()};
    std::set<unsigned int> bounce_back_dirs = bounce_back::determine_bounce_back_directions(current_border_info);
    
    for (const auto i : bounce_back_dirs)
    {
        remaining_dirs.erase(i);
    }
    unsigned int x = std::get<0>(lbm_access::get_node_coordinates(current_border_info[0]));
    unsigned int y = std::get<1>(lbm_access::get_node_coordinates(current_border_info[0]));
    if(x == 1)
    {
            // if(y == 1) remaining_dirs.insert({2,5});
            // else if(y == (VERTICAL_NODES - 2)) remaining_dirs.insert({5,8});
            // else remaining_dirs.insert({2,5,8});
        remaining_dirs.insert({2,5,8});
    }
    else if(x ==(HORIZONTAL_NODES - 2))
    {
        // if(y == 1) remaining_dirs.insert({0,3});
        // else if(y == (VERTICAL_NODES - 2)) remaining_dirs.insert({3,6});
        // else remaining_dirs.insert({0,3,6});
        remaining_dirs.insert({0,3,6});
    }
    return remaining_dirs;  
}
