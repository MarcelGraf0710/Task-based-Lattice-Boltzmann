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
namespace bounce_back
{
    
    /**
     * @brief Retrieves the border adjacencies for all fluid nodes within the simulation domain based on 
     *        the phase information of all nodes. Notice that all fluid nodes on the edges of the simulation 
     *        domain will automatically become border nodes.
     * 
     * @param fluid_nodes a vector containing the indices of all fluid nodes within the simulation domain
     * @param phase_information a vector containing the phase information of ALL nodes 
     * @return the border adjancency relations of each node as it is required by bounce-back boundary treatment
     *   
     */
    border_adjacency retrieve_border_adjacencies
    (
        std::vector<unsigned int> &fluid_nodes, 
        std::vector<bool> &phase_information
    );

    /**
     * @brief Retrieves the border swap information for all fluid nodes within the simulation domain based on the 
     *        phase information of all nodes. Notice that all fluid nodes on the edges of the simulation domain will 
     *        automatically become border nodes.
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
     * @brief Determines the directions in which the value will be reflected from the opposite direction.
     * 
     * @param current_border_info an entry of a border_swap_information object
     * @return a set containing all bounce back directions
     */
    std::set<unsigned int> determine_bounce_back_directions
    (
        std::vector<unsigned int> &current_border_info
    );

    /**
     * @brief Performs a halfway bounce-back streaming update for all fluid nodes within the simulation domain.
     *        This version utilizes the ghost nodes bordering a boundary node. It is intended for use with
     *        the two-step, swap and shift algorithms.
     * 
     * @param ba see documentation of border_adjacency
     * @param distribution_values a vector containing the distribution values of all nodes
     * @param access_function the access function used to access the distribution values
     */
    void perform_boundary_update
    (
        border_adjacency &ba,
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
    void perform_early_boundary_update
    (
        border_swap_information &border_nodes,
        std::vector<double> &source, // NEW!!
        std::vector<double> &destination, // RENAMED
        access_function access_function
    );

    /**
     * @brief Modified version of the halfway bounce-back streaming update for all fluid nodes 
     *        within the simulation domain. Instead of using information stored in ghost nodes, 
     *        This allows for a convenient unification of the streaming and collision step for
     *        the two-lattice algorithm.
     *        This variant will update inlets and outlets according to the specified velocity and density.
     * 
     * @param bsi see documentation of border_swap_information
     * @param distribution_values a vector containing the distribution values of all nodes
     * @param velocities a vector containing all velocities
     * @param densities a vector containing all densities
     * @param access_function the access function used to access the distribution values
     */
    void perform_inout_boundary_update
    (
        border_swap_information &bsi,
        std::vector<double> &distribution_values, 
        // std::vector<velocity> &velocities,
        // std::vector<double> densities,
        access_function access_function
    );

    /**
     * @brief Updates the ghost nodes that represent inlet and outlet edges
     * 
     * @param distribution_values a vector containing the distribution values of all nodes
     * @param access_function the access function used to access the distribution values
     */
    void inout_update
    (
        std::vector<double> &distribution_values, 
        access_function access_function
    );
}

#endif