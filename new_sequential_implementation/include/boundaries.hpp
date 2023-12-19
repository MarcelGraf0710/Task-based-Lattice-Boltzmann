#ifndef BOUNDARIES_HPP
#define BOUNDARIES_HPP
#include "defines.hpp"

/**
 * @brief Returns whether the node located at the edge of the simulation domain.
 * 
 * @param node_index the index of the node in question
 * @return true if the coordinates are any of the node are (1, y), (HORIZONTAL_NODES - 2, y), (x, 1), (x, VERTICAL_NODES - 2)
 *         with suitable x and y
 */
bool is_edge_node(unsigned int node_index);

/**
 * @brief Returns whether the node is a ghost node
 * 
 * @param node_index the index of the node in question
 */
bool is_ghost_node(unsigned int node_index);

/**
 * @brief Returns whether the node is located at the edge of the domain
 * 
 * @param x x coordinate of the node in question
 * @param y y coordinate of the node in question
 * @return true if the coordinates are any of the node are (1, y), (HORIZONTAL_NODES - 2, y), (x, 1), (x, VERTICAL_NODES - 2)
 *         with suitable x and y
 */
inline bool is_edge_node(unsigned int x, unsigned int y)
{
    bool result = ((x == 1) || (x = HORIZONTAL_NODES - 2)) && ((y == 1) || (y = VERTICAL_NODES - 2));
    return result;
}

/**
 * @brief Returns a vector containing all fluid non-border nodes within the simulation domain
 * 
 * @param fluid_nodes a vector containing all fluid nodes within the simulation domain
 * @param ba see documentation of border_adjacency
 */
std::vector<unsigned int> get_non_border_nodes
(
    std::vector<unsigned int> &fluid_nodes,
    border_adjacency ba
);

/**
 * @brief This namespace contains all function representations of boundary conditions used in the lattice-Boltzmann model.
 */
namespace bounce_back
{
    
    /**
     * @brief Retrieves the border adjacencies for all fluid nodes within the simulation domain based on the phase information of all nodes.
     *        Notice that all fluid nodes on the edges of the simulation domain will automatically become border nodes.
     * 
     * @param fluid_nodes a vector containing the indices of all fluid nodes within the simulation domain
     * @param phase_information a vector containing the phase information of ALL nodes 
     * @return the border adjancency relations of each node as it is required by bounce-back boundary treatment
     */
    border_adjacency retrieve_border_adjacencies
    (
        std::vector<unsigned int> &fluid_nodes, 
        std::vector<bool> &phase_information
    );

    /**
     * @brief Retrieves the border swap information for all fluid nodes within the simulation domain based on the phase information of all nodes.
     *        Notice that all fluid nodes on the edges of the simulation domain will automatically become border nodes.
     * 
     * @param fluid_nodes a vector containing the indices of all fluid nodes within the simulation domain
     * @param phase_information a vector containing the phase information of ALL nodes 
     * @return see documentation of border_swap_information
     */
    border_swap_information retrieve_border_swap_information
    (
        std::vector<unsigned int> &fluid_nodes, 
        std::vector<bool> &phase_information
    );

    /**
     * @brief Performs a bounce-back streaming update for all fluid nodes within the simulation domain.
     * 
     * @param border_nodes see documentation of border_adjacency
     * @param distribution_values a vector containing the distribution values of all nodes
     * @param access_function the access function used to access the distribution values
     */
    void perform_boundary_update(
        border_adjacency &border_nodes,
        std::vector<double> &distribution_values, 
        access_function access_function
        );

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
    void perform_early_boundary_update(
        border_swap_information &border_nodes,
        std::vector<double> &distribution_values, 
        access_function access_function
    );

}

#endif