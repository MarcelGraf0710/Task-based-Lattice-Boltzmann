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
typedef std::vector<std::tuple<unsigned int, char, char>> swap_information;

namespace swap_sequential
{
    /**
     * @brief Maps the border node types to the swap partner directions.
     *        0 - {}        (upper right corner node)
     *        1 - {5}       (upper edge node or left upper corner node)
     *        2 - {6,7}     (right edge node or right lower corner node)
     *        3 - {5,7,8}   (left edge node or left lower corner node)
     *        4 - {5,6,7,8} (non-border node or lower edge node)
     */
    extern std::map<char, std::set<unsigned int>> swap_directions;

    /**
     * @brief This map stores the inlet and outlet directions for each node.
     *        0 - {}        (empty set, neither inlet nor outlet)
     *        1 - {2,5,8}   (inlet node)
     *        2 - {0,3,6}   (outlet node)
     */
    extern std::map<char, std::set<unsigned int>> inout_directions;

    /**
     * @brief Sets up the swap information data structure
     * 
     * @param bsi see documentation of border_swap_information
     * @param fluid_nodes a vector containing the indices of all fluid nodes
     * @return see documentation of swap_information 
     */
    swap_information setup_swap_information
    (
        const border_swap_information &bsi, 
        const std::vector<unsigned int> &fluid_nodes
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
        unsigned int node_index,
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
     * @param swap_directions a set containing all directions in which swap will take place
     */
    inline void perform_swap_step
    (
        std::vector<double> &distribution_values,
        unsigned int node_index,
        const access_function access_function,
        const std::set<unsigned int> &swap_directions
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
     * @brief Performs the sequential swap algorithm for the specified number of iterations.
     * 
     * @param fluid_nodes A vector containing the indices of all fluid nodes in the domain
     * @param values the vector containing the distribution values of all nodes
     * @param access_function the access function according to which the values are to be accessed
     * @param iterations this many iterations will be performed
     */
    void run
    (  
        std::vector<unsigned int> &fluid_nodes,       
        std::vector<double> &values, 
        border_swap_information &bsi,
        access_function access_function,
        unsigned int iterations
    );

    /**
     * @brief Performs the combined streaming and collision step for all fluid nodes within the simulation domain.
     *        The border conditions are enforced through ghost nodes.
     *        This variant of the combined streaming and collision step will print several debug comments to the console.
     */
    sim_data_tuple perform_swap_stream_and_collide_debug
    (
        const swap_information &swap_information,
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
        const swap_information &swap_information,
        std::vector<double> &distribution_values,    
        const access_function access_function
    );
}