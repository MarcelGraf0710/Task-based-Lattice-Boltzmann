#include "../include/swap_sequential.hpp"
#include "../include/access.hpp"
#include "../include/boundaries.hpp"
#include "../include/macroscopic.hpp"
#include "../include/collision.hpp"
#include <iostream>

const std::vector<unsigned int> swap_sequential::ACTIVE_STREAMING_DIRECTIONS = {5,6,7,8};

/**
 * @brief Retrieves an improved version of the border swap information data structure.
 *        This method does not consider inlet and outlet ghost nodes when performing bounce-back
 *        as the inserted values will be overwritten by inflow and outflow values anyways.
 * 
 * @param fluid_nodes a vector containing the indices of all fluid nodes within the simulation domain
 * @param phase_information a vector containing the phase information for every vector (true means solid)
 * @return border_swap_information see documentation of border_swap_information
 */
border_swap_information swap_sequential::retrieve_swap_info
(
    const std::vector<unsigned int> &fluid_nodes, 
    const std::vector<bool> &phase_information
)
{
    border_swap_information result;
    std::vector<unsigned int> current_adjacencies;
    std::vector<unsigned int> swap_adjacencies;
    std::tuple<unsigned int, unsigned int> coords;
    for(const auto node : fluid_nodes)
    {
        current_adjacencies = {};
        swap_adjacencies = {};
        // Determine all nodes with non-inout ghost neighbors
        for(const auto direction : STREAMING_DIRECTIONS)
        {
            unsigned int current_neighbor = access::get_neighbor(node, direction);
            if(is_non_inout_ghost_node(current_neighbor, phase_information))
            {
                current_adjacencies.push_back(direction);
            }
        }
        // Determine intersection with active streaming direction mask
        std::set_intersection(
            current_adjacencies.begin(), 
            current_adjacencies.end(), 
            ACTIVE_STREAMING_DIRECTIONS.begin(), 
            ACTIVE_STREAMING_DIRECTIONS.end(), 
            std::back_inserter(swap_adjacencies)
            );

        // Determine node coordinates
        coords = access::get_node_coordinates(node);
        
        if(std::get<0>(coords) == 1) // Inlet node
        {
            swap_adjacencies.push_back(0);
            swap_adjacencies.push_back(3);
        }
        else if(std::get<0>(coords) == (HORIZONTAL_NODES - 2)) // Outlet node
        {
            swap_adjacencies.push_back(2);
        }
        current_adjacencies = {node};
        std::sort(swap_adjacencies.begin(), swap_adjacencies.end());
        current_adjacencies.insert(current_adjacencies.end(), swap_adjacencies.begin(), swap_adjacencies.end());
        
        if(current_adjacencies.size() > 1) result.push_back(current_adjacencies);
    }
    return result;
}

/**
 * @brief Performs the combined streaming and collision step for all fluid nodes within the simulation domain.
 *        The border conditions are enforced through ghost nodes.
 *        This variant of the combined streaming and collision step will print several debug comments to the console.
 */
sim_data_tuple swap_sequential::perform_swap_stream_and_collide_debug
(
    const border_swap_information &bsi,
    const std::vector<unsigned int> &fluid_nodes,
    std::vector<double> &distribution_values,    
    const access_function access_function
)
{
    std::cout << "Distribution values before stream and collide: " << std::endl;
    to_console::print_distribution_values(distribution_values, access_function);
    std::cout << std::endl;

    std::vector<velocity> velocities(TOTAL_NODE_COUNT, velocity{0,0});
    velocity current_velocity = {0,0};

    std::vector<double> densities(TOTAL_NODE_COUNT, -1);
    double current_density = -1;

    unsigned int node_index = 0;
    std::vector<unsigned int> swap_dirs;
    std::vector<double> current_distributions(DIRECTION_COUNT, 0);

    // Border node initialization
    for(const auto node : bsi)
    {
        for(auto it = node.begin() + 1; it < node.end(); ++it)
        {
            swap_sequential::perform_swap_step(distribution_values, node[0], access_function, *it);
        }
    }
    std::cout << "Distribution values after border node initialization: " << std::endl;
    to_console::print_distribution_values(distribution_values, access_function);

    // Swapping step
    for(const auto node : fluid_nodes)
    {
        swap_sequential::perform_swap_step(distribution_values, node, access_function, ACTIVE_STREAMING_DIRECTIONS);
    }
    std::cout << "Distribution values after swap for every node: " << std::endl;
    to_console::print_distribution_values(distribution_values, access_function);

    // Restore precious order in here
    for(const auto node : fluid_nodes)
    {
        swap_sequential::restore_order(distribution_values, node, access_function);
    }
    std::cout << "Distribution values after ORDER has been restored for every node: " << std::endl;
    to_console::print_distribution_values(distribution_values, access_function);


    // Collision
    for(const auto& node : fluid_nodes)
    {
        current_distributions = access::get_distribution_values_of(distribution_values, node, access_function);
        current_velocity = macroscopic::flow_velocity(current_distributions);
        current_density = macroscopic::density(current_distributions);
        velocities[node] = current_velocity;
        densities[node] = current_density;
        current_distributions = collision::collide_bgk(current_distributions, current_velocity, current_density);
        access::set_distribution_values_of(current_distributions, distribution_values, node, access_function);
    }
    std::cout << "Distribution values after collision: " << std::endl;
    to_console::print_distribution_values(distribution_values, access_function);

    /* Update ghost nodes */
    boundary_conditions::update_velocity_input_density_output(distribution_values, velocities, densities, access_function);
    std::cout << "Distribution values after ghost node update: " << std::endl;
    to_console::print_distribution_values(distribution_values, access_function);
    std::cout << std::endl;

    sim_data_tuple result{velocities, densities};

    return result;
}

/**
 * @brief Performs the combined streaming and collision step for all fluid nodes within the simulation domain.
 *        The border conditions are enforced through ghost nodes.
 *        This variant of the combined streaming and collision step will print several debug comments to the console.
 */
sim_data_tuple swap_sequential::perform_swap_stream_and_collide
(
    const border_swap_information &bsi,
    const std::vector<unsigned int> &fluid_nodes,
    std::vector<double> &distribution_values,    
    const access_function access_function
)
{
    std::vector<velocity> velocities(TOTAL_NODE_COUNT, velocity{0,0});
    velocity current_velocity = {0,0};

    std::vector<double> densities(TOTAL_NODE_COUNT, -1);
    double current_density = -1;

    std::vector<double> current_distributions(DIRECTION_COUNT, 0);

    // Border node initialization
    for(const auto node : bsi)
    {
        for(auto it = node.begin() + 1; it < node.end(); ++it)
        {
            swap_sequential::perform_swap_step(distribution_values, node[0], access_function, *it);
        }
    }

    for(const auto node : fluid_nodes)
    {
        // Swapping step
        swap_sequential::perform_swap_step(distribution_values, node, access_function, ACTIVE_STREAMING_DIRECTIONS);

        // Restore precious order in here
        swap_sequential::restore_order(distribution_values, node, access_function);

        // Collision
        current_distributions = access::get_distribution_values_of(distribution_values, node, access_function);
        current_velocity = macroscopic::flow_velocity(current_distributions);
        current_density = macroscopic::density(current_distributions);
        velocities[node] = current_velocity;
        densities[node] = current_density;
        current_distributions = collision::collide_bgk(current_distributions, current_velocity, current_density);
        access::set_distribution_values_of(current_distributions, distribution_values, node, access_function);
    }

    /* Update ghost nodes */
    boundary_conditions::update_velocity_input_density_output(distribution_values, velocities, densities, access_function);

    sim_data_tuple result{velocities, densities};

    return result;
}

/**
 * @brief Performs the sequential swap algorithm for the specified number of iterations.
 * 
 * @param fluid_nodes A vector containing the indices of all fluid nodes in the domain
 * @param values the vector containing the distribution values of all nodes
 * @param access_function the access function according to which the values are to be accessed
 * @param iterations this many iterations will be performed
 */
void swap_sequential::run
(  
    const std::vector<unsigned int> &fluid_nodes,
    const std::vector<bool> &phase_information,       
    std::vector<double> &values, 
    const access_function access_function,
    const unsigned int iterations
)
{
    to_console::print_run_greeting("sequential swap algorithm", iterations);    

    std::vector<sim_data_tuple>result(
        iterations, 
        std::make_tuple(std::vector<velocity>(TOTAL_NODE_COUNT, {0,0}), std::vector<double>(TOTAL_NODE_COUNT, 0)));

    std::cout << "Retrieving swap info" << std::endl;
    border_swap_information bsi = swap_sequential::retrieve_swap_info(fluid_nodes, phase_information);

    for(auto time = 0; time < iterations; ++time)
    {
        std::cout << "\033[33mIteration " << time << ":\033[0m" << std::endl;
        result[time] = swap_sequential::perform_swap_stream_and_collide(bsi, fluid_nodes, values, access_function);
        std::cout << "\tFinished iteration " << time << std::endl;
    }

    to_console::print_simulation_results(result);
    std::cout << "All done, exiting simulation. " << std::endl;
}