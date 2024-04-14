#include "../include/two_lattice_parallel.hpp"
#include "../include/access.hpp"
#include "../include/defines.hpp"
#include "../include/boundaries.hpp"
#include "../include/collision.hpp"
#include "../include/macroscopic.hpp"
#include "../include/utils.hpp"
#include "../include/parallel_framework.hpp"
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
    parallel_framework::emplace_bounce_back_values(bsi, source, access_function);

    /* Combined stream and collision step */
    hpx::for_each
    (
        hpx::execution::par, 
        fluid_nodes.begin(), 
        fluid_nodes.end(), 
        [&source, &destination, access_function, &velocities, &densities](unsigned int fluid_node)
        {
            tl_stream_and_collide_helper(source, destination, access_function, fluid_node, velocities, densities);
        }
    );

    two_lattice_parallel::update_velocity_input_density_output(destination, velocities, densities, access_function);

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
    two_lattice_parallel::tl_stream(
        source,
        destination,
        access_function,
        fluid_node);

    std::vector<double> current_distributions =
        lbm_access::get_distribution_values_of(destination, fluid_node, access_function);

    velocity current_velocity = macroscopic::flow_velocity(current_distributions);
    velocities[fluid_node] = current_velocity;

    double current_density = macroscopic::density(current_distributions);
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
    std::vector<double> densities(TOTAL_NODE_COUNT, -1);

    std::cout << "\t TL stream and collide: initializations and declarations performed." << std::endl;

    /* Boundary node treatment */
    bounce_back::emplace_bounce_back_values(bsi, source, access_function);
    std::cout << "SOURCE after emplace bounce-back values: " << std::endl;
    to_console::print_distribution_values(source, access_function);

    /* Streaming step */
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

    std::cout << "DESTINATION after streaming: " << std::endl;
    to_console::print_distribution_values(destination, access_function);

    /* Collision step */
    hpx::for_each
    (
        hpx::execution::par, 
        fluid_nodes.begin(), 
        fluid_nodes.end(), 
        [&velocities, &densities, &source, &destination, access_function](unsigned int fluid_node)
        {
            std::vector<double> current_distributions = 
                lbm_access::get_distribution_values_of(destination, fluid_node, access_function);

            velocity current_velocity = macroscopic::flow_velocity(current_distributions);    
            velocities[fluid_node] = current_velocity;

            double current_density = macroscopic::density(current_distributions);
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
        result[time] = two_lattice_parallel::perform_tl_stream_and_collide
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
 * @brief Updates the ghost nodes that represent inlet and outlet edges.
 *        When updating, a velocity border condition will be considered for the input
 *        and a density border condition for the output.
 *        The inlet velocity is constant throughout all inlet nodes whereas the outlet nodes
 *        all have the specified density.
 *        The corresponding values are constants defined in "../include/"defines.hpp".
 * 
 * @param distribution_values a vector containing the distribution values of all nodes
 * @param velocities a vector containing the velocities of all nodes
 * @param densities a vector containing the densities of all nodes
 * @param access_function the access function used to access the distribution values
 */
void two_lattice_parallel::update_velocity_input_density_output
(
    std::vector<double> &distribution_values,
    std::vector<velocity> &velocities,
    std::vector<double> &densities, 
    const access_function access_function
)
{
    hpx::experimental::for_loop(
        hpx::execution::par, 0, VERTICAL_NODES - 1,
        [&distribution_values, &velocities, &densities, access_function](int y)
        {
        // Update inlets
        int current_border_node = lbm_access::get_node_index(0,y);
        velocity v = INLET_VELOCITY;
        double density = INLET_DENSITY;
        std::vector<double> current_dist_vals = maxwell_boltzmann_distribution(v, density);
        lbm_access::set_distribution_values_of
        (
            current_dist_vals,
            distribution_values,
            current_border_node,
            access_function
        );
        velocities[current_border_node] = v;
        densities[current_border_node] = density;

        // Update outlets
        current_border_node = lbm_access::get_node_index(HORIZONTAL_NODES - 1,y);
        v = macroscopic::flow_velocity(lbm_access::get_distribution_values_of(distribution_values, lbm_access::get_neighbor(current_border_node, 3), access_function));
        density = OUTLET_DENSITY;
        current_dist_vals = maxwell_boltzmann_distribution(v, density);
        lbm_access::set_distribution_values_of
        (
            current_dist_vals,
            distribution_values,
            current_border_node,
            access_function
        );
        velocities[current_border_node] = v;
        densities[current_border_node] = density;
        });
}