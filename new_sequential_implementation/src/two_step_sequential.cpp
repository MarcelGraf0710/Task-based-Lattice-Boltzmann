#include "../include/two_step_sequential.hpp"
#include "../include/access.hpp"
#include "../include/two_step_sequential.hpp"
#include "../include/boundaries.hpp"
#include "../include/macroscopic.hpp"
#include "../include/new_collision.hpp"

/**
 * @brief Performs a vertical streaming step for all fluid nodes, i.e. all f_i with i in {1,7} will be propagated.
 * 
 * @param fluid_nodes A vector containing the indices of all fluid nodes in the domain
 * @param values a vector containing all distribution values.
 * @param access_function the access to node values will be performed according to this access function.
 * @param positive_direction true if the velocity vector index is greater than four, false otherwise
 */
void two_step_sequential::perform_vertical_stream
(
    std::vector<unsigned int> &fluid_nodes,       
    std::vector<double> &values, 
    access_function access_function
)
{
    // Downward streaming with direction 1
    for(auto it = fluid_nodes.begin() + 1; it < fluid_nodes.end() - 2; ++it)
    {
        values[access_function(access::get_neighbor(*it, 1), 1)] = values[access_function(*it, 1)];
    }

    // Upward streaming with direction 7
    for(auto it = fluid_nodes.end() - 2; it >= fluid_nodes.begin(); --it)
    {
        values[access_function(access::get_neighbor(*it, 7), 7)] = values[access_function(*it, 7)];
    }
}

/**
 * @brief Performs a horizontal streaming step for all fluid nodes, i.e. all f_i with i in {3,5} will be propagated.
 * 
 * @param fluid_nodes A vector containing the indices of all fluid nodes in the domain
 * @param values a vector containing all distribution values.
 * @param access_function the access to node values will be performed according to this access function.
 * @param positive_direction true if the velocity vector index is greater than four, false otherwise
 */
void two_step_sequential::perform_horizontal_stream
(
    std::vector<unsigned int> &fluid_nodes,       
    std::vector<double> &values, 
    access_function access_function
)
{
    unsigned int current_index = 0;

    // Leftward streaming with direction 3
    for(auto column = 1; column < HORIZONTAL_NODES-1; ++column)
    {
        for(auto row = 1; row < VERTICAL_NODES-1; ++row)
        {
            current_index = access::get_node_index(column, row);
            values[access_function(access::get_neighbor(current_index, 3), 3)] = values[access_function(current_index, 3)];
        }
    }

    // Rightward streaming with direction 5
    for(auto column = HORIZONTAL_NODES-1; column >= 1 ; --column)
    {
        for(auto row = 1; row < VERTICAL_NODES-1; ++row)
        {
            current_index = access::get_node_index(column, row);
            values[access_function(access::get_neighbor(current_index, 5), 5)] = values[access_function(current_index, 5)];
        }
    }
}

/**
 * @brief Performs the streaming step for all fluid nodes within the simulation domain.
 * 
 * @param fluid_non_border_nodes A vector containing the indices of all fluid non-border nodes in the domain
 * @param ba see documentation of border_adjacency
 * @param values a vector containing all distribution values
 * @param access_function the access to node values will be performed according to this access function.
 */
void two_step_sequential::perform_fast_stream
(
    std::vector<unsigned int> &fluid_non_border_nodes, 
    border_adjacency ba,
    std::vector<double> &values, 
    access_function access_function
)
{
    // All directions that require left-to-right and/or bottom-to-top node iteration order
    for(auto it = fluid_non_border_nodes.begin() + 1; it < fluid_non_border_nodes.end() - 1; ++it)
    {
        values[access_function(access::get_neighbor(*it, 0), 0)] = values[access_function(*it, 0)];
        values[access_function(access::get_neighbor(*it, 1), 1)] = values[access_function(*it, 1)];
        values[access_function(access::get_neighbor(*it, 2), 2)] = values[access_function(*it, 2)];
        values[access_function(access::get_neighbor(*it, 3), 3)] = values[access_function(*it, 3)];
    }

    // All directions that require right-to-left and/or top-to-bottom node iteration order
    for(auto it = fluid_non_border_nodes.end() - 2; it > fluid_non_border_nodes.begin(); --it)
    {
        values[access_function(access::get_neighbor(*it, 5), 5)] = values[access_function(*it, 5)];
        values[access_function(access::get_neighbor(*it, 6), 6)] = values[access_function(*it, 6)];
        values[access_function(access::get_neighbor(*it, 7), 7)] = values[access_function(*it, 7)];
        values[access_function(access::get_neighbor(*it, 8), 8)] = values[access_function(*it, 8)];
    }

    // Update boundary values
    bounce_back::perform_boundary_update(ba, values, access_function);
}

/**
 * @brief Performs the sequential two-step algorithm for the specified number of iterations.
 * 
 * @param fluid_nodes A vector containing the indices of all fluid nodes in the domain
 * @param values the vector containing the distribution values of all nodes
 * @param access_function the access function according to which the values are to be accessed
 * @param iterations this many iterations will be performed
 */
std::vector<sim_data_tuple> run
(  
    std::vector<unsigned int> &fluid_nodes,       
    std::vector<double> &values, 
    border_adjacency &ba,
    access_function access_function,
    unsigned int iterations
)
{
    std::cout << "Now running sequential two lattice algorithm for " << iterations << " iterations." << std::endl;
    // Create vector that contains only non-border nodes
    // For the speficied number of iterations...
    // Perform the streaming step for all nodes
    // Calculate all macroscopic sizes
    // Perform the collision step for all nodes
    // Repeat
    std::vector<unsigned int> non_border_nodes = get_non_border_nodes(fluid_nodes, ba);
    std::vector<sim_data_tuple> sim_data;
    std::vector<double> current_density = {};
    for(auto i = 0; i < iterations; ++i)
    {
        std::cout << "\t Iteration " << time << ":" << std::endl;
        two_step_sequential::perform_fast_stream(non_border_nodes, ba, values, access_function);
        std::cout << "\t\t Streaming step performed " << std::endl;
        sim_data.push_back(macroscopic::get_sim_data_tuple(fluid_nodes, values, access_function));
        current_density = std::get<1>(sim_data[i]);
        collision::collide_all_bgk(fluid_nodes, values, std::get<0>(sim_data[i]), std::get<1>(sim_data[i]), access_function);
        std::cout << "\t\t Collision step performed " << std::endl;
    }
}
