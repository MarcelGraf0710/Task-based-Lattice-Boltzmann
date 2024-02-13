#ifndef BOUNDARIES_HPP
#define BOUNDARIES_HPP
#include "defines.hpp"
#include <set>

/**
 * @brief Returns whether the node with the specified index is located at the edge of the simulation domain.
 *        This is the case for any of the following coordinates:
 *        - (1, y)
 *        - (HORIZONTAL_NODES - 2, y)
 *        - (x, 1)
 *        - (x, VERTICAL_NODES - 2)
 *        with suitable x and y.
 * 
 * @param node_index the index of the node in question
 * @return whether or not the node is an edge node
 *         
 */
bool is_edge_node(unsigned int node_index);

/**
 * @brief Returns whether the node with the specified index is a ghost node.
 *        This is the case for any of the following:
 *        
 *        With suitable x and y, the node coordinates are any of
 *        - (0, y)
 *        - (HORIZONTAL_NODES - 1, y)
 *        - (x, 0)
 *        - (x, VERTICAL_NODES - 1)
 *        
 *        or
 * 
 *        The node with the specified index is solid.
 * 
 * @param node_index the index of the node in question
 * @param phase_information a vector containing the phase information for all nodes of the lattice
 */
bool is_ghost_node(unsigned int node_index, std::vector<bool> &phase_information);

/**
 * @brief Returns a vector containing all fluid non-border nodes within the simulation domain.
 * 
 * @param fluid_nodes a vector containing all fluid nodes within the simulation domain
 * @param ba see documentation of border_adjacency
 */
std::vector<unsigned int> get_non_border_nodes
(
    std::vector<unsigned int> &fluid_nodes,
    border_adjacency &ba
);

/**
 * @brief This namespace contains all function representations of boundary conditions used in the lattice-Boltzmann model.
 */
namespace velocity_boundary
{
    void update_lower_wall_boundaries
    (
        std::vector<double> &source,
        std::vector<double> &destination, 
        access_function &access_function, 
        std::set<unsigned int> &directions
    )
    {
        // for(auto direction : directions) 
        // {
        //     //std::cout << "\t performing stream (node = " << fluid_node << ", dir = " << direction << ") := " << "(node = " << access::get_neighbor(fluid_node, invert_direction(direction)) << ", dir = " << direction << ")" << std::endl; 
        //     destination[access_function(fluid_node, direction)] =
        //     source[
        //         access_function(
        //             access::get_neighbor(fluid_node, invert_direction(direction)), 
        //             direction)];
        // }
    }
}

namespace outlet
{

}
#endif