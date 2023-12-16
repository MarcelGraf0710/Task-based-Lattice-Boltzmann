#include "../include/new_two_lattice.hpp"
#include "../include/access.hpp"
#include "../include/defines.hpp"
#include "../include/boundaries.hpp"
#include "../include/new_collision.hpp"
#include "../include/macroscopic.hpp"
#include <set>

/**
 * @brief Performs the combined streaming-and-collision step for all nodes within the simulation domain.
 * 
 * @param time_step current iteration of the simulation
 * @param fluid_nodes A vector containing the indices of all fluid nodes in the domain
 * @param boundary_nodes A vector containing the indices of all fluid boundary nodes in the domain
 * @param values_0 source for even time steps and destination for odd time steps
 * @param values_1 source for odd time steps and destination for even time steps
 * @param access_function the access function according to which the values are to be accessed
 */
sim_data_tuple two_lattice_sequential::perform_tl_stream_and_collide
(
    std::vector<unsigned int> &fluid_nodes,
    border_swap_information &boundary_nodes,
    std::vector<double> &source, 
    std::vector<double> &destination,    
    access_function access_function
)
{
    std::set<unsigned int> remaining_nodes = {fluid_nodes.begin(), fluid_nodes.end()};
    std::set<unsigned int> remaining_dirs;
    std::vector<double> vals;
    velocity current_velocity;
    double density;
    std::vector<velocity> velocities(fluid_nodes.size(), velocity{0,0});
    std::vector<double> densities(fluid_nodes.size(), 0);

    // Pre-streaming boundary condition fulfillment
    bounce_back::perform_early_boundary_update(boundary_nodes, destination, access_function);

    // For every boundary node, update remaining streaming directions
    for(auto current : boundary_nodes)
    {
        // Determine remaining directions
        remaining_dirs = {streaming_directions.begin(), streaming_directions.end()};
        for (auto i = current.begin() + 1; i< current.end(); ++i)
        {
            remaining_dirs.erase(*i);
        }
        for(auto dir : remaining_dirs)
        {
            destination[access_function(current[0], dir)] = 
            source[access_function(access::get_neighbor(current[0], dir), invert_direction(dir))];
        }

        // Perform collision
        vals = access::get_all_distribution_values(destination, current[0], access_function);
        current_velocity = macroscopic::flow_velocity(vals);
        velocities[current[0]] = current_velocity;
        density = macroscopic::density(vals);
        densities[current[0]] = density;
        vals = collision::collide_bgk(vals, current_velocity, density);
        access::set_all_distribution_values(vals, destination, current[0], access_function);
        
        // Delete ready boundary node from remaining nodes
        remaining_nodes.erase(current[0]);
    }

    // Perform streaming and collision for remaining fluid nodes
    for(auto fluid_node : remaining_nodes)
    {
        // Perform streaming
        for (auto direction = 0; direction < DIRECTION_COUNT; ++direction)
        {
            destination[access_function(fluid_node, direction)] = 
            source[access_function(access::get_neighbor(fluid_node, direction), invert_direction(direction))];

        }

        // Perform collision
        vals = access::get_all_distribution_values(destination, fluid_node, access_function);
        current_velocity = macroscopic::flow_velocity(vals);
        velocities[fluid_node] = current_velocity;
        density = macroscopic::density(vals);
        densities[fluid_node] = density;
        vals = collision::collide_bgk(vals, current_velocity, density);
        access::set_all_distribution_values(vals, destination, fluid_node, access_function);
    }

    return sim_data_tuple{velocities, densities};
}

/**
 * @brief Performs the sequential two-lattice algorithm for the specified number of iterations.
 * 
 * @param fluid_nodes A vector containing the indices of all fluid nodes in the domain
 * @param boundary_nodes A vector containing the indices of all fluid boundary nodes in the domain
 * @param values_0 source for even time steps and destination for odd time steps
 * @param values_1 source for odd time steps and destination for even time steps
 * @param access_function the access function according to which the values are to be accessed
 * @param iterations this many iterations will be performed
 * @return a vector of tuples containing all flow velocities and density values for all time steps
 */
std::vector<sim_data_tuple> two_lattice_sequential::run
(  
    std::vector<unsigned int> &fluid_nodes,       
    border_swap_information &boundary_nodes,
    std::vector<double> &values_0, 
    std::vector<double> &values_1,   
    access_function access_function,
    unsigned int iterations
)
{
    std::vector<double> &source = values_0;
    std::vector<double> &destination = values_1;
    std::vector<double> &temp = values_1;
    std::vector<sim_data_tuple> result;

    for(auto time = 0; time < iterations; ++time)
    {
        result.push_back(
            two_lattice_sequential::perform_tl_stream_and_collide(
                fluid_nodes, boundary_nodes, source, destination, access_function)
            );
        
        temp = source;
        source = destination;
        destination = temp;
    }
    
    return result;
}


