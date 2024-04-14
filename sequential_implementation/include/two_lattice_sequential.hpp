#ifndef TWO_LATTICE_SEQUENTIAL_HPP
#define TWO_LATTICE_SEQUENTIAL_HPP

#include <vector>
#include <set>
#include "defines.hpp"
#include "access.hpp"
#include <iostream>
#include "collision.hpp"
#include "utils.hpp"
#include "boundaries.hpp"

namespace sequential_two_lattice
{
    /**
     * @brief Performs the combined streaming and collision step for all fluid nodes within the simulation domain.
     *        The border conditions are enforced through ghost nodes.
     * 
     * @param fluid_nodes a vector containing the indices of all fluid nodes within the simulation domain.
     * @param bsi see documentation of border_swap_information
     * @param source a vector containing the distribution values of the previous time step
     * @param destination the distribution values will be written to this vector after performing both steps.
     * @param access_function the function used to access the distribution values
     * @return see documentation of sim_data_tuple
     */
    sim_data_tuple perform_tl_stream_and_collide
    (
        const std::vector<unsigned int> &fluid_nodes,
        const border_swap_information &bsi,
        std::vector<double> &source, 
        std::vector<double> &destination,    
        const access_function access_function
    );

    /**
     * @brief Performs the combined streaming and collision step for all fluid nodes within the simulation domain.
     *        The border conditions are enforced through ghost nodes.
     *        This variant of the combined streaming and collision step will print several debug comments to the console.
     * 
     * @param fluid_nodes a vector containing the indices of all fluid nodes within the simulation domain.
     * @param bsi see documentation of border_swap_information
     * @param source a vector containing the distribution values of the previous time step
     * @param destination the distribution values will be written to this vector after performing both steps.
     * @param access_function the function used to access the distribution values
     * @return see documentation of sim_data_tuple
     */
    sim_data_tuple perform_tl_stream_and_collide_debug
    (
        const std::vector<unsigned int> &fluid_nodes,
        const border_swap_information &bsi,
        std::vector<double> &source, 
        std::vector<double> &destination,    
        const access_function access_function
    );

    /**
     * @brief Performs the sequential two-lattice algorithm for the specified number of iterations.
     * 
     * @param fluid_nodes A vector containing the indices of all fluid nodes in the domain
     * @param boundary_nodes see documentation of border_swap_information
     * @param distribution_values_0 source for even time steps and destination for odd time steps
     * @param distribution_values_1 source for odd time steps and destination for even time steps
     * @param access_function the access function according to which distribution values are to be accessed
     * @param iterations this many iterations will be performed
     */
    void run
    (  
        const std::vector<unsigned int> &fluid_nodes,       
        const border_swap_information &boundary_nodes,
        std::vector<double> &distribution_values_0, 
        std::vector<double> &distribution_values_1,   
        const access_function access_function,
        const unsigned int iterations
    );

    /**
     * @brief Determines the remaining streaming option for a node based on the specified border 
     *        information vector.
     * 
     * @param current_border_info an entry of a border_swap_information object
     * @return a set containing all remaining streaming directions
     */
    std::set<unsigned int> determine_streaming_directions
    (
        const std::vector<unsigned int> &current_border_info
    );

    /**
     * @brief Performs the steaming step in all directions for the fluid node with 
     *        the specified index.
     * 
     * @param source distribution values will be taken from this vector
     * @param destination distribution values will be rearranged in this vector
     * @param access_function function that will be used to access the distribution values
     * @param fluid_node the index of the node for which the streaming step is performed
     */
    inline void tl_stream
    (
        const std::vector<double> &source,
        std::vector<double> &destination, 
        const access_function &access_function, 
        const unsigned int fluid_node
    )
    {
        for (const auto direction : ALL_DIRECTIONS)
        {
            destination[access_function(fluid_node, direction)] =
                source[
                    access_function(
                        lbm_access::get_neighbor(fluid_node, invert_direction(direction)), 
                        direction)];
        }
    }

    /**
     * @brief Performs the collision step for the fluid node with the specified index.
     * 
     * @param destination the updated distribution values will be written to this vector
     * @param fluid_node the index of the fluid node
     * @param access_function the function used to access the distribution values
     * @param velocity the velocity at the node in question
     * @param density the density at the node in question
     */
    inline void tl_collision
    (
        std::vector<double> &destination, 
        const unsigned int fluid_node, 
        const std::vector<double> &distribution_values,
        const access_function &access_function, 
        const velocity &velocity, 
        const double &density
    )
    {
        std::vector<double> vals = collision::collide_bgk(distribution_values, velocity, density);
        lbm_access::set_distribution_values_of(vals, destination, fluid_node, access_function);
    }
}

#endif