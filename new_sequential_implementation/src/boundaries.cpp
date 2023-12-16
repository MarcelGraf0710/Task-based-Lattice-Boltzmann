#include "../include/boundaries.hpp"
#include "../include/access.hpp"
#include "../include/defines.hpp"
#include <iostream>

/**
 * @brief Returns whether the node located at the edge of the simulation domain.
 * 
 * @param node_index the index of the node in question
 * @return true if the coordinates are any of the node are (1, y), (HORIZONTAL_NODES - 2, y), (x, 1), (x, VERTICAL_NODES - 2)
 *         with suitable x and y
 */
bool is_edge_node(unsigned int node_index)
{
    std::tuple<unsigned int, unsigned int> coordinates = access::get_node_coordinates(node_index);
    unsigned int x = std::get<0>(coordinates);
    unsigned int y = std::get<1>(coordinates);
    bool result = ((x == 1) || (x = HORIZONTAL_NODES - 2)) && ((y == 1) || (y = VERTICAL_NODES - 2));
    return result;
}

/**
 * @brief Returns whether the node is a ghost node
 * 
 * @param node_index the index of the node in question
 */
bool is_ghost_node(unsigned int node_index)
{
    std::tuple<unsigned int, unsigned int> coordinates = access::get_node_coordinates(node_index);
    unsigned int x = std::get<0>(coordinates);
    unsigned int y = std::get<1>(coordinates);
    bool result = ((x == 0) || (x = HORIZONTAL_NODES - 1)) && ((y == 0) || (y = VERTICAL_NODES - 1));
    return result;
}

/**
 * @brief Retrieves the border adjacencies for all fluid nodes within the simulation domain based on the phase information of all nodes.
 *        Notice that all fluid nodes on the edges of the simulation domain will automatically become border nodes.
 * 
 * @param fluid_nodes a vector containing the indices of all fluid nodes within the simulation domain
 * @param phase_information a vector containing the phase information of ALL nodes 
 * @return the border adjancency relations of each node as it is required by bounce-back boundary treatment
 */
border_adjacency bounce_back::retrieve_border_adjacencies
(        
    std::vector<unsigned int> &fluid_nodes, 
    std::vector<bool> &phase_information
)
{
    border_adjacency result;

    for(auto node : fluid_nodes)
    {
        std::vector<std::tuple<unsigned int, unsigned int>> current_adjacencies{std::make_tuple(node, 4)};
        for(auto direction : {0,1,2,3,5,6,7,8})
        {
            unsigned int current_neighbor = access::get_neighbor(node, direction);
            if(is_ghost_node(current_neighbor) | phase_information[current_neighbor])
            {
                current_adjacencies.push_back(std::make_tuple(current_neighbor, direction));
            }
            if(current_adjacencies.size() > 1) result.push_back(current_adjacencies);
        }
    }
    return result;
}

/**
 * @brief Retrieves the border swap information for all fluid nodes within the simulation domain based on the phase information of all nodes.
 *        Notice that all fluid nodes on the edges of the simulation domain will automatically become border nodes.
 * 
 * @param fluid_nodes a vector containing the indices of all fluid nodes within the simulation domain
 * @param phase_information a vector containing the phase information of ALL nodes 
 * @return see documentation of border_swap_information
 */
border_swap_information bounce_back::retrieve_border_swap_information
(
    std::vector<unsigned int> &fluid_nodes, 
    std::vector<bool> &phase_information
)
{
    border_swap_information result;

    for(auto node : fluid_nodes)
    {
        std::vector<unsigned int> current_adjacencies{node};
        for(auto direction : {0,1,2,3,5,6,7,8})
        {
            unsigned int current_neighbor = access::get_neighbor(node, direction);
            if(is_ghost_node(current_neighbor) | phase_information[current_neighbor])
            {
                current_adjacencies.push_back(direction);
            }
            if(current_adjacencies.size() > 1) result.push_back(current_adjacencies);
        }
    }
    return result;
}

/**
 * @brief Performs a bounce-back streaming update for all fluid nodes within the simulation domain.
 * 
 * @param border_nodes see documentation of border_adjacency
 * @param distribution_values a vector containing the distribution values of all nodes
 * @param access_function the access function used to access the distribution values
 */
void bounce_back::perform_boundary_update
(
    border_adjacency &border_nodes,
    std::vector<double> &distribution_values, 
    access_function access_function
)
{
    for(auto current_border_node_information : border_nodes)
    {
        int current_border_node = std::get<0>(current_border_node_information[0]);
        std::vector<std::tuple<unsigned int, unsigned int>> neighbors
        {
            current_border_node_information.begin() + 1, 
            current_border_node_information.end()
        };
        for(auto neighbor : neighbors)
        {
            unsigned int current_dir = invert_direction(std::get<1>(neighbor));
            distribution_values[access_function(current_border_node, current_dir)] = 
            distribution_values[access_function(std::get<0>(neighbor), current_dir)];
        }
    }
}

/**
 * @brief Modified version of the bounce-back streaming update for all fluid nodes within the simulation domain.
 *        Instead of using information stored in ghost nodes, it inverts the pre-streaming values such that it is 
 *        already available when a streaming step concluded. This allows for the collision step to be merged with
 *        the streaming step.
 * 
 * @param border_nodes see documentation of border_swap_information
 * @param distribution_values a vector containing the distribution values of all nodes
 * @param access_function the access function used to access the distribution values
 */
void bounce_back::perform_early_boundary_update
(
    border_swap_information &border_nodes,
    std::vector<double> &distribution_values, 
    access_function access_function
)
{
    for(auto current : border_nodes)
    {
        int current_border_node = current[0];
        for(auto direction = current.begin() + 1; direction < current.end(); ++direction)
        {
            distribution_values[access_function(current_border_node, invert_direction(*direction))] = 
            distribution_values[access_function(current_border_node, current[*direction])];
        }
    }
}
