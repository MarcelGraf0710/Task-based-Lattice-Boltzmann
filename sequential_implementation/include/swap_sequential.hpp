#ifndef TWO_SWAP_SEQUENTIAL_HPP
#define TWO_SWAP_SEQUENTIAL_HPP

#include <map>
#include <set>
#include <vector>
#include "defines.hpp"
#include "utils.hpp"

/**
 * @brief For each fluid node, this data structure stores 
 *        - the index of the node
 *        - the key for the swap directions 
 *        - the key for the inout directions
 */
typedef std::vector<std::tuple<unsigned int, std::set<unsigned int>>> swap_information;

namespace swap_sequential
{
    extern const std::vector<unsigned int> ACTIVE_STREAMING_DIRECTIONS;

    border_swap_information retrieve_swap_info
    (
        const std::vector<unsigned int> &fluid_nodes, 
        const std::vector<bool> &phase_information
    );

    /**
     * @brief Restores the original order after streaming is completed for a node.
     *        Notice that this method automatically achieves the halfway bounce-back boundary condition.
     * 
     * @param distribution_values a vector containing the distribution values of all nodes
     * @param node_index the index of the node whose values are to be reordered
     * @param access_function the function used to access the distribution values
     */
    inline void restore_order
    (
        std::vector<double> &distribution_values,
        const unsigned int node_index,
        const access_function access_function
    )
    {
        for(const auto dir : {0,1,2,3})
        {
            vec_utils::swap(
                distribution_values, 
                access_function(node_index, dir), 
                access_function(node_index, invert_direction(dir))
            );
        }
    }

    /**
     * @brief Performs the swap step that is equvalent to the "active" streaming step for the specified node.
     * 
     * @param distribution_values a vector containing the distribution values of all nodes
     * @param node_index the index of the node for which the streaming step is to be performed
     * @param access_function the function used to access the distribution values
     * @param swap_directions a vector containing all directions in which swap will take place
     */
    inline void perform_swap_step
    (
        std::vector<double> &distribution_values,
        const unsigned int node_index,
        const access_function access_function,
        const std::vector<unsigned int> &swap_directions
    )
    {
        for(const auto dir : swap_directions)
        {
            vec_utils::swap(
                distribution_values, 
                access_function(node_index, dir), 
                access_function(access::get_neighbor(node_index, dir), invert_direction(dir))
            );
        }
    }

    /**
     * @brief Performs the swap step in a single direction for the specified node.
     * 
     * @param distribution_values a vector containing the distribution values of all nodes
     * @param node_index the index of the node for which the streaming step is to be performed
     * @param access_function the function used to access the distribution values
     * @param direction the direction in which swap will take place
     */
    inline void perform_swap_step
    (
        std::vector<double> &distribution_values,
        const unsigned int node_index,
        const access_function access_function,
        const unsigned int direction
    )
    {
        vec_utils::swap(
            distribution_values, 
            access_function(node_index, direction), 
            access_function(access::get_neighbor(node_index, direction), invert_direction(direction))
        );
    }

    /**
     * @brief Performs the sequential swap algorithm for the specified number of iterations.
     * 
     * @param fluid_nodes A vector containing the indices of all fluid nodes in the domain
     * @param values the vector containing the distribution values of all nodes
     * @param access_function the access function according to which the values are to be accessed
     * @param iterations this many iterations will be performed
     */
    void run
    (  
        const std::vector<unsigned int> &fluid_nodes,
        const std::vector<bool> &phase_information,       
        std::vector<double> &values, 
        const access_function access_function,
        const unsigned int iterations
    );

    /**
     * @brief Performs the combined streaming and collision step for all fluid nodes within the simulation domain.
     *        The border conditions are enforced through ghost nodes.
     *        This variant of the combined streaming and collision step will print several debug comments to the console.
     */
    sim_data_tuple perform_swap_stream_and_collide_debug
    (
        const border_swap_information &bsi,
        const std::vector<unsigned int> &fluid_nodes,
        std::vector<double> &distribution_values,    
        const access_function access_function
    );

    /**
     * @brief Performs the combined streaming and collision step for all fluid nodes within the simulation domain.
     *        The border conditions are enforced through ghost nodes.
     *        This variant of the combined streaming and collision step will print several debug comments to the console.
     */
    sim_data_tuple perform_swap_stream_and_collide
    (
        const border_swap_information &bsi,
        const std::vector<unsigned int> &fluid_nodes,
        std::vector<double> &distribution_values,    
        const access_function access_function
    );
}

#endif