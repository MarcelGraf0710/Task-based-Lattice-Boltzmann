#ifndef BOUNDARIES_HPP
#define BOUNDARIES_HPP
#include "defines.hpp"
#include "access.hpp"
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
bool is_ghost_node(unsigned int node_index, const std::vector<bool> &phase_information);

/**
 * @brief Returns a vector containing all fluid non-border nodes within the simulation domain.
 * 
 * @param fluid_nodes a vector containing all fluid nodes within the simulation domain
 * @param ba see documentation of border_adjacency
 */
std::vector<unsigned int> get_non_border_nodes
(
    const std::vector<unsigned int> &fluid_nodes,
    const border_adjacency &ba
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
        const std::vector<unsigned int> &fluid_nodes, 
        const std::vector<bool> &phase_information
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
        const std::vector<unsigned int> &fluid_nodes, 
        const std::vector<bool> &phase_information
    );

    /**
     * @brief Determines the directions in which the value will be reflected from the opposite direction.
     * 
     * @param current_border_info an entry of a border_swap_information object
     * @return a set containing all bounce back directions
     */
    std::set<unsigned int> determine_bounce_back_directions
    (
        const std::vector<unsigned int> &current_border_info
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
        const border_swap_information &bsi,
        std::vector<double> &distribution_values, 
        const access_function access_function
    );

    /**
     * @brief Modified version of the halfway bounce-back streaming update for all fluid nodes 
     *        within the simulation domain. Instead of using information stored in ghost nodes, 
     *        This allows for a convenient unification of the streaming and collision step for
     *        the two-lattice algorithm.
     * 
     * @param bsi see documentation of border_swap_information
     * @param source the distribution values will be read from this vector
     * @param destination the updated distribution values will be written to this vector
     * @param access_function the access function used to access the distribution values
     */
    void perform_early_boundary_update
    (
        const border_swap_information &bsi,
        const std::vector<double> &source, 
        std::vector<double> &destination, 
        const access_function access_function
    );
}

/**
 * @brief This namespace contains all functions that are required for enforcing boundary conditions.
 *        It uses the ghost nodes for this purpose.
 * 
 */
namespace boundary_conditions
{
    /**
     * @brief Modified version of the halfway bounce-back streaming update for all fluid nodes 
     *        within the simulation domain. Instead of using information stored in ghost nodes, 
     *        This allows for a convenient unification of the streaming and collision step for
     *        the two-lattice algorithm.
     *        This variant will update inlets and outlets according to the specified velocity and density.
     * 
     * @param bsi see documentation of border_swap_information
     * @param distribution_values a vector containing the distribution values of all nodes
     * @param access_function the access function used to access the distribution values
     */
    void perform_inout_boundary_update
    (
        const border_swap_information &bsi,
        std::vector<double> &distribution_values, 
        const access_function access_function
    );

    /**
     * @brief Updates the ghost nodes that represent inlet and outlet edges
     * 
     * @param distribution_values a vector containing the distribution values of all nodes
     * @param access_function the access function used to access the distribution values
     */
    void update_velocity_input_velocity_output
    (
        std::vector<double> &distribution_values, 
        const access_function access_function
    );

    /**
     * @brief Updates the ghost nodes that represent inlet and outlet edges.
     *        When updating, a velocity border condition will be considered for the input
     *        and a density border condition for the output.
     * 
     * @param distribution_values a vector containing the distribution values of all nodes
     * @param access_function the access function used to access the distribution values
     */
    void update_velocity_input_density_output
    (
        std::vector<double> &distribution_values, 
        const access_function access_function
    );

    /**
     * @brief Updates the ghost nodes that represent inlet and outlet edges.
     *        When updating, a density border condition will be considered for both the input and the output.
     *        The corresponding values are constants defined in "../include/"defines.hpp".
     * 
     * @param distribution_values a vector containing the distribution values of all nodes
     * @param access_function the access function used to access the distribution values
     */
    void update_density_input_density_output
    (
        std::vector<double> &distribution_values, 
        const access_function access_function
    );

    /**
     * @brief Initializes all inlet and outlet nodes with their corresponding initial values.
     *        The corresponding values are constants defined in "../include/"defines.hpp".
     * 
     * @param distribution_values a vector containing the distribution values of all nodes
     * @param access_function the access function used to access the distribution values
     */
    void initialize_inout
    (
        std::vector<double> &distribution_values, 
        const access_function access_function
    );

    /**
     * @brief Realizes inflow and outflow by an inward stream of each border node.
     *        This method is intended for use with two-step, swap and shift algorithms.
     * 
     * @param distribution_values a vector containing the distribution values of all nodes
     * @param access_function the access function used to access the distribution values
     */
    void ghost_stream_inout
    (
        std::vector<double> &distribution_values, 
        const access_function access_function
    );

    /**
     * @brief Performs inflow or outflow for the specified node in the specified directions
     * 
     * @param distribution_values a vector containing all distribution values
     * @param node_index the index of the node for which inflow or outflow is to be performed
     * @param dirs {2,5,8} if node_index is an inlet 
     *             {0,3,6} if node_index is an outlet
     *             {} if node_index is neither an inlet nor an outlet
     * @param access_function 
     */
    inline void single_node_inout
    (
    std::vector<double> &distribution_values,
    unsigned int node_index,
    const std::set<unsigned int> &dirs,
    const access_function access_function
    )
    {
        for(const auto direction : dirs)
        {
            distribution_values[access_function(node_index, direction)] = 
                distribution_values[access_function(access::get_neighbor(node_index, invert_direction(direction)), direction)];
        }
    }
}

/**
 * @brief This namespace contains discrete versions of important velocity profiles that occur for streaming fluids.
 * 
 */
namespace velocity_profiles
{
    /**
     * @brief Computes a laminary velocity profile for inlet or outlet nodes.
     * 
     * @param u the mean velocity of the profile
     * @return std::vector<velocity> a vector containing the velocity values for the inlet or outlet nodes.
     */
    std::vector<velocity> ideal_laminary(velocity &u);

    /**
     * @brief Computes a turbulent velocity profile for inlet or outlet nodes using the rule of the seventh.
     * 
     * @param u the mean velocity of the profile
     * @return std::vector<velocity> a vector containing the velocity values for the inlet or outlet nodes.
     */
    std::vector<velocity> seventh_rule_turbulent(velocity &u);
}

#endif