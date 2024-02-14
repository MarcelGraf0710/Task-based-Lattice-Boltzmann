#include "../include/new_two_lattice.hpp"
#include "../include/access.hpp"
#include "../include/defines.hpp"
#include "../include/boundaries.hpp"
#include "../include/new_collision.hpp"
#include "../include/macroscopic.hpp"
#include "../include/utils.hpp"
#include <set>
#include <iostream>


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
sim_data_tuple two_lattice_sequential::perform_tl_stream_and_collide
(
    const std::vector<unsigned int> &fluid_nodes,
    const border_swap_information &bsi,
    const std::vector<double> &source, 
    std::vector<double> &destination,    
    const access_function access_function
)
{
    std::set<unsigned int> remaining_nodes = {fluid_nodes.begin(), fluid_nodes.end()};
    std::set<unsigned int> remaining_dirs = {streaming_directions.begin(), streaming_directions.end()};
    std::vector<velocity> velocities(TOTAL_NODE_COUNT, velocity{0,0});
    std::vector<double> densities(TOTAL_NODE_COUNT, -1);
    std::vector<double> current_distributions(DIRECTION_COUNT, 0);
    std::vector<double> current;

    /* Boundary node treatment */
    bounce_back::perform_early_boundary_update(bsi, source, destination, access_function);

    /* Combined streaming and collision step */
    for(auto current_border_info : bsi)
    {
        remaining_dirs = determine_streaming_directions(current_border_info);
        two_lattice_sequential::tl_stream(
            source, 
            destination, 
            access_function, 
            current_border_info[0], 
            remaining_dirs);

        current_distributions = 
            access::get_distribution_values_of(destination, current_border_info[0], access_function);
        velocities[current_border_info[0]] = macroscopic::flow_velocity(current_distributions);
        densities[current_border_info[0]] = macroscopic::density(current_distributions);

        two_lattice_sequential::tl_collision(
            destination, 
            current_border_info[0], 
            current_distributions,
            access_function, 
            velocities, 
            densities);

        remaining_nodes.erase(current_border_info[0]);
    }

    /* Treatment of non-boundary nodes */
    for(auto fluid_node : remaining_nodes)
    {   
        two_lattice_sequential::tl_stream(
            source, 
            destination, 
            access_function, 
            fluid_node);

        current_distributions = 
            access::get_distribution_values_of(destination, fluid_node, access_function);
            
        velocities[fluid_node] = macroscopic::flow_velocity(current_distributions);
        densities[fluid_node] = macroscopic::density(current_distributions);
        
        two_lattice_sequential::tl_collision(
            destination, 
            fluid_node, 
            current_distributions,
            access_function, 
            velocities, 
            densities);
    }

    /* Update ghost nodes */
    boundary_conditions::update_velocity_input_density_output(destination, access_function);

    unsigned int update_node = 0;
    for(auto x = 0; x < HORIZONTAL_NODES; x = x + HORIZONTAL_NODES - 1)
    {
        for(auto y = 1; y < VERTICAL_NODES - 1; ++y)
        {
            update_node = access::get_node_index(x,y);
            current_distributions = 
                access::get_distribution_values_of(destination, update_node, access_function);
                
            velocities[update_node] = macroscopic::flow_velocity(current_distributions);
            densities[update_node] = macroscopic::density(current_distributions);
        }
    }

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
sim_data_tuple two_lattice_sequential::perform_tl_stream_and_collide_debug
(
    const std::vector<unsigned int> &fluid_nodes,
    const border_swap_information &bsi,
    const std::vector<double> &source, 
    std::vector<double> &destination,    
    const access_function access_function
)
{
    std::cout << "\t SOURCE before stream and collide: " << std::endl;
    to_console::print_distribution_values(source, access_function);
    std::cout << std::endl;

    std::cout << "DESTINATION as received by perform_tl_stream_and_collide: " << std::endl;
    to_console::print_distribution_values(destination, access_function);
    std::cout << std::endl;
    destination = source;

    std::set<unsigned int> remaining_nodes = {fluid_nodes.begin(), fluid_nodes.end()};
    std::set<unsigned int> remaining_dirs = {streaming_directions.begin(), streaming_directions.end()};
    std::vector<velocity> velocities(TOTAL_NODE_COUNT, velocity{0,0});
    std::vector<double> densities(TOTAL_NODE_COUNT, -1);
    std::vector<double> current_distributions(DIRECTION_COUNT, 0);
    std::vector<double> current;

    std::cout << "\t TL stream and collide: initializations and declarations performed." << std::endl;

    /* Boundary node treatment */

    bounce_back::perform_early_boundary_update(bsi, source, destination, access_function);
    std::cout << "\t Early boundary update performed." << std::endl;

    std::cout << "SOURCE after boundary update: " << std::endl;
    to_console::print_distribution_values(source, access_function);
    std::cout << std::endl;

    std::cout << "DESTINATION after boundary update: " << std::endl;
    to_console::print_distribution_values(destination, access_function);
    std::cout << std::endl;

    for(auto current_border_info : bsi)
    {
        remaining_dirs = determine_streaming_directions(current_border_info);
        two_lattice_sequential::tl_stream(
            source, 
            destination, 
            access_function, 
            current_border_info[0], 
            remaining_dirs);
    }

    std::cout << "DESTINATION after combined stream: " << std::endl;
    to_console::print_distribution_values(destination, access_function);
    std::cout << std::endl;

    for(auto current_border_info : bsi)
    {
        current_distributions = 
            access::get_distribution_values_of(destination, current_border_info[0], access_function);
        velocities[current_border_info[0]] = macroscopic::flow_velocity(current_distributions);
        densities[current_border_info[0]] = macroscopic::density(current_distributions);

        two_lattice_sequential::tl_collision(
            destination, 
            current_border_info[0], 
            current_distributions,
            access_function, 
            velocities, 
            densities);

        remaining_nodes.erase(current_border_info[0]);
    }

    std::cout << "DESTINATION after combined collide: " << std::endl;
    to_console::print_distribution_values(destination, access_function);
    std::cout << std::endl;

    std::cout << "\t Performed stream and collision for all border nodes." << std::endl;

    /* Treatment of non-boundary nodes */
    for(auto fluid_node : remaining_nodes)
    {   
        two_lattice_sequential::tl_stream(
            source, 
            destination, 
            access_function, 
            fluid_node);
    }

    for(auto fluid_node : remaining_nodes)
    {   
        current_distributions = 
            access::get_distribution_values_of(destination, fluid_node, access_function);
            
        velocities[fluid_node] = macroscopic::flow_velocity(current_distributions);
        densities[fluid_node] = macroscopic::density(current_distributions);
        
        two_lattice_sequential::tl_collision(
            destination, 
            fluid_node, 
            current_distributions,
            access_function, 
            velocities, 
            densities);
    }
    std::cout << "\t Done treating all non-boundary nodes." << std::endl;
    std::cout << "\t Distribution values after streaming and collision: " << std::endl;
    to_console::print_distribution_values(destination, access_function);
    std::cout << std::endl;
    std::cout << "\t PASSED #1 " << std::endl;

    boundary_conditions::update_density_input_density_output(destination, access_function);
    std::cout << "Updated inlet and outlet ghost nodes... hopefully. Check it out:" <<std::endl;
    to_console::print_distribution_values(destination, access_function);

    unsigned int update_node = 0;
    for(auto x = 0; x < HORIZONTAL_NODES; x = x + HORIZONTAL_NODES - 1)
    {
        for(auto y = 1; y < VERTICAL_NODES - 1; ++y)
        {
            update_node = access::get_node_index(x,y);
            current_distributions = 
                access::get_distribution_values_of(destination, update_node, access_function);
                
            velocities[update_node] = macroscopic::flow_velocity(current_distributions);
            densities[update_node] = macroscopic::density(current_distributions);
        }
    }

    sim_data_tuple result{velocities, densities};

    return result;
}

/**
 * @brief Performs the sequential two-lattice algorithm for the specified number of iterations.
 * 
 * @param fluid_nodes A vector containing the indices of all fluid nodes in the domain
 * @param boundary_nodes A vector containing the indices of all fluid boundary nodes in the domain
 * @param values_0 source for even time steps and destination for odd time steps
 * @param values_1 source for odd time steps and destination for even time steps
 * @param access_function the access function according to which distribution values are to be accessed
 * @param iterations this many iterations will be performed
 * @param data the simulation data tuples will be placed in this vector (assumed to be pre-initialized)
 */
void two_lattice_sequential::run
(  
    const std::vector<unsigned int> &fluid_nodes,       
    const border_swap_information &boundary_nodes,
    std::vector<double> &values_0, 
    std::vector<double> &values_1,   
    const access_function access_function,
    unsigned int iterations,
    std::vector<sim_data_tuple> &data
)
{
    std::cout << "------------------------------------------------------------------------------------------------------------------------" << std::endl;
    std::cout << "Now running sequential two lattice algorithm for " << iterations << " iterations." << std::endl;
    std::cout << std::endl;

    std::vector<double> &source = values_0;
    std::vector<double> &destination = values_1;
    std::vector<double> &temp = values_1;
    std::vector<sim_data_tuple>result(
        iterations, 
        std::make_tuple(std::vector<velocity>(TOTAL_NODE_COUNT, {0,0}), std::vector<double>(TOTAL_NODE_COUNT, 0)));

    for(auto time = 0; time < iterations; ++time)
    {
        std::cout << "\033[33mIteration " << time << ":\033[0m" << std::endl;
        result[time] = two_lattice_sequential::perform_tl_stream_and_collide
        (
            fluid_nodes, 
            boundary_nodes, 
            source, 
            destination, 
            access_function
        );
            
        std::cout << "Finished iteration " << time << std::endl;
        std::swap(source, destination);
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

/**
 * @brief Determines the remaining streaming option for a node based on the specified border 
 *        information vector.
 * 
 * @param current_border_info an entry of a border_swap_information object
 * @return a set containing all remaining streaming directions
 */
std::set<unsigned int> two_lattice_sequential::determine_streaming_directions
(
    const std::vector<unsigned int> &current_border_info
)
{
    std::set<unsigned int> remaining_dirs = {streaming_directions.begin(), streaming_directions.end()};
    std::set<unsigned int> bounce_back_dirs = bounce_back::determine_bounce_back_directions(current_border_info);
    
    for (auto i : bounce_back_dirs)
    {
        remaining_dirs.erase(i);
    }
    unsigned int x = std::get<0>(access::get_node_coordinates(current_border_info[0]));
    unsigned int y = std::get<1>(access::get_node_coordinates(current_border_info[0]));
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
