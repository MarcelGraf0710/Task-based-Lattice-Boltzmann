#include "../include/new_two_lattice.hpp"
#include "../include/access.hpp"
#include "../include/defines.hpp"
#include "../include/boundaries.hpp"

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
void two_lattice_sequential::perform_tl_stream_and_collide
    (
        std::vector<double> &fluid_nodes,
        border_swap_information &boundary_nodes,
        std::vector<double> &source, 
        std::vector<double> &destination,    
        access_function access_function
    )
    {
        // Pre-streaming boundary condition fulfillment
        bounce_back::perform_early_boundary_update(boundary_nodes, destination, access_function);

        for(auto fluid_node : fluid_nodes)
        {
            // Streaming executed for all fluid nodes, including boundary nodes
            for (auto direction = 0; direction < DIRECTION_COUNT; ++direction)
            {
                destination[access_function(access::get_neighbor(fluid_node, direction), direction)] = source[access_function(fluid_node, direction)];
            }
            
        }
        
        


    }
