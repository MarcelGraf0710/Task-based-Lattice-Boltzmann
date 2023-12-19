#include "../include/new_two_lattice.hpp"
#include "../include/access.hpp"
#include "../include/defines.hpp"
#include "../include/boundaries.hpp"
#include "../include/new_collision.hpp"
#include "../include/macroscopic.hpp"
#include <set>
#include <stdexcept>

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
    sim_data_tuple result{velocities, densities};

    // Pre-streaming boundary condition fulfillment
    bounce_back::perform_early_boundary_update(boundary_nodes, destination, access_function);
    std::cout << "\t Early boundary update performed." << std::endl;

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

    std::cout << "\t Performed stream and collision for all border nodes." << std::endl;

    // Perform streaming and collision for remaining fluid nodes
    for(auto fluid_node : remaining_nodes)
    {
        //std::cout << "Currently dealing with node " << fluid_node;
        /* TEMP */ std::tuple<int, int> coords = access::get_node_coordinates(fluid_node);
        //std::cout << ", coords: (" << std::get<0>(coords) << ", " << std::get<1>(coords) << ")" << std::endl; 

        // Perform streaming
        for (auto direction = 0; direction < DIRECTION_COUNT; ++direction)
        {
            destination[access_function(fluid_node, direction)] = 
            source[access_function(access::get_neighbor(fluid_node, direction), invert_direction(direction))];

        }
        //std::cout << "Streaming performed. " << std::endl;

        // Perform collision
        vals = access::get_all_distribution_values(destination, fluid_node, access_function);
        //if(fluid_node == 681) std::cout << "Determined all distribution values " << std::endl;
        current_velocity = macroscopic::flow_velocity(vals);
        //if(fluid_node == 681) std::cout << "Determined current velocity " << std::endl;
        velocities[fluid_node] = current_velocity;
        density = macroscopic::density(vals);
        //if(fluid_node == 681) std::cout << "Determined density " << std::endl;
        densities[fluid_node] = density;
        //if(fluid_node == 681) std::cout << "About to access collide_bgk " << std::endl;
        vals = collision::collide_bgk(vals, current_velocity, density);
        //if(fluid_node == 681) std::cout << "Collision bgk accessed " << std::endl;
        access::set_all_distribution_values(vals, destination, fluid_node, access_function);

        //std::cout << "Collision performed. " << std::endl;
    }
    std::cout << "\t Done." << std::endl;
    //std::cout << std::get<1>(result)[30] << std::endl;
    //return sim_data_tuple{velocities, densities};
    return result;
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
    std::cout << "Now running sequential two lattice algorithm for " << iterations << " iterations." << std::endl;
    std::vector<double> &source = values_0;
    std::vector<double> &destination = values_1;
    std::vector<double> &temp = values_1;
    std::vector<sim_data_tuple> result;

    for(auto time = 0; time < iterations; ++time)
    {
        std::cout << "Iteration " << time << ":" << std::endl;
        result.push_back(
            two_lattice_sequential::perform_tl_stream_and_collide(
                fluid_nodes, boundary_nodes, source, destination, access_function)
            );
        std::cout << "\tBoth stream and collide performed for all nodes, back at run(...), now changing source and destination..." << std::endl;
        temp = source;
        source = destination;
        destination = temp;
    }
    std::cout << "All done. " << std::endl;
    
    return result;
}


