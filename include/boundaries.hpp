#ifndef BOUNDARIES_HPP
#define BOUNDARIES_HPP

#include "access.hpp"
#include "defines.hpp"

#include <set>

/**
 * @brief This set contains all directions that point inside the simulation domain for a regular inlet node.
 */
extern std::set<unsigned int> INFLOW_INSTREAM_DIRS;

/**
 * @brief This set contains all directions that point inside the simulation domain for a regular outlet node.
 */
extern std::set<unsigned int> OUTFLOW_INSTREAM_DIRS;

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
inline bool is_edge_node(unsigned int node_index)
{
    std::tuple<unsigned int, unsigned int> coordinates = lbm_access::get_node_coordinates(node_index);
    unsigned int x = std::get<0>(coordinates);
    unsigned int y = std::get<1>(coordinates);
    bool result = ((x == 1) || (x == (HORIZONTAL_NODES - 2))) && ((y == 1) || (y == (VERTICAL_NODES - 2)));
    return result;
}

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
inline bool is_ghost_node(unsigned int node_index, const std::vector<bool> &phase_information)
{
    std::tuple<unsigned int, unsigned int> coordinates = lbm_access::get_node_coordinates(node_index);
    unsigned int x = std::get<0>(coordinates);
    unsigned int y = std::get<1>(coordinates);
    bool is_outer_node = ((x == 0) || (x == (HORIZONTAL_NODES - 1))) || ((y == 0) || (y == (VERTICAL_NODES - 1)));
    return is_outer_node  | phase_information[node_index];
}

/**
 * @brief Determines whether a node with the specified index is a non-inlet and non-outlet ghost node.
 *        This holds true if all of the following are fulfilled:
 *        - The node is not located at the absolute left or right of the lattice
 *        - The node is located at the absolute top or bottom of the lattice OR it is solid
 * 
 * @param node_index the index of the node in question
 * @param phase_information a vector storing the phase information of all lattice nodes
 * @return true if the node is neither an inlet node nor an outlet node, and false otherwise
 */
inline bool is_non_inout_ghost_node
(
    unsigned int node_index, 
    const std::vector<bool> &phase_information
)
{
    std::tuple<unsigned int, unsigned int> coordinates = lbm_access::get_node_coordinates(node_index);
    unsigned int x = std::get<0>(coordinates);
    unsigned int y = std::get<1>(coordinates);
    bool is_outer_non_inout_node = ((x != 0) && (x != (HORIZONTAL_NODES - 1))) && ((y == 0) || (y == (VERTICAL_NODES - 1)) || phase_information[node_index]);
    return is_outer_non_inout_node;
}

/**
 * @brief This namespace contains function representations of boundary conditions used in the lattice-Boltzmann model.
 *        Notice that parallel versions exist within "parallel_framework.hpp".
 */
namespace bounce_back
{

    /**
     * @brief Retrieves the border swap information data structure.
     *        This method does not consider inlet and outlet ghost nodes when performing bounce-back
     *        as the inserted values will be overwritten by inflow and outflow values anyways.
     * 
     * @param fluid_nodes a vector containing the indices of all fluid nodes within the simulation domain
     * @param phase_information a vector containing the phase information for every vector (true means solid)
     * @return border_swap_information see documentation of border_swap_information
     */
    border_swap_information retrieve_border_swap_info
    (
        const std::vector<unsigned int> &fluid_nodes, 
        const std::vector<bool> &phase_information
    );

    /**
     * @brief Performs an outstream step for all border nodes in the directions where they border non-inout ghost nodes.
     *        The distribution values will be stored in the ghost nodes in inverted order such that
     *        after this method is executed, the border nodes can be treated like regular nodes when performing an instream.
     * 
     * @param bsi a border_swap_information generated by retrieve_fast_border_swap_info
     * @param distribution_values a vector containing the distribution values of all nodes
     * @param access_function the access function used to access the distribution values
     * @param read_offset offset for shift algorithm, leave zero for all other algorithms
     */
    void emplace_bounce_back_values
    (
        const border_swap_information &bsi,
        std::vector<double> &distribution_values,
        const access_function access_function,
        const unsigned int read_offset = 0
    );

    /**
     * @brief Performs a halfway bounce-back streaming update for all fluid nodes within the simulation domain.
     *        This version utilizes the ghost nodes bordering a boundary node. It is intended for use with
     *        the two-step algorithm.
     * 
     * @param bsi see documentation of border_swap_information
     * @param distribution_values a vector containing the distribution values of all nodes
     * @param access_function the access function used to access the distribution values
     */
    void perform_boundary_update
    (
        const border_swap_information &bsi,
        std::vector<double> &distribution_values, 
        const access_function access_function
    );
}

/**
 * @brief This namespace contains all functions that are required for enforcing boundary conditions.
 *        It uses the ghost nodes for this purpose.
 */
namespace boundary_conditions
{
    /**
     * @brief Updates the ghost nodes that represent inlet and outlet edges.
     *        When updating, a velocity border condition will be considered for both the input and the output.
     *        The corresponding values are constants defined in "defines.hpp".
     * 
     * @param distribution_values a vector containing the distribution values of all nodes
     * @param velocities a vector containing the velocities of all nodes
     * @param densities a vector containing the densities of all nodes
     * @param access_function the access function used to access the distribution values
     */
    void update_velocity_input_velocity_output
    (
        std::vector<double> &distribution_values,
        std::vector<velocity> &velocities,
        std::vector<double> &densities, 
        const access_function access_function
    );

    /**
     * @brief Updates the ghost nodes that represent inlet and outlet edges.
     *        When updating, a velocity border condition will be considered for the input
     *        and a density border condition for the output.
     *        The inlet velocity is constant throughout all inlet nodes whereas the outlet nodes
     *        all have the specified density.
     *        The corresponding values are constants defined in "../include/"defines.hpp".
     * 
     * @param distribution_values a vector containing the distribution values of all nodes
     * @param velocities a vector containing the velocities of all nodes
     * @param densities a vector containing the densities of all nodes
     * @param access_function the access function used to access the distribution values
     */
    void update_velocity_input_density_output
    (
    std::vector<double> &distribution_values,
    std::vector<velocity> &velocities,
    std::vector<double> &densities, 
    const access_function access_function
    );

    /**
     * @brief Updates the ghost nodes that represent inlet and outlet edges.
     *        When updating, a density border condition will be considered for both the input and the output.
     *        The corresponding values are constants defined in "../include/"defines.hpp".
     * 
     * @param distribution_values a vector containing the distribution values of all nodes
     * @param velocities a vector containing the velocities of all nodes
     * @param densities a vector containing the densities of all nodes
     * @param access_function the access function used to access the distribution values
     */
    void update_density_input_density_output
    (
        std::vector<double> &distribution_values, 
        std::vector<velocity> &velocities,
        std::vector<double> &densities, 
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
     *        This method is intended for use with the two-step algorithm.
     * 
     * @param distribution_values a vector containing the distribution values of all nodes
     * @param access_function the access function used to access the distribution values
     */
    void ghost_stream_inout
    (
        std::vector<double> &distribution_values, 
        const access_function access_function
    );
}

/**
 * @brief This namespace contains discrete versions of important velocity profiles that occur for streaming fluids.
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