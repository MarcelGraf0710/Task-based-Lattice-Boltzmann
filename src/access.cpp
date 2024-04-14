#include "../include/access.hpp"
#include <list>

/**
* @brief Retrieves the coordinates of the node with the specified node index.
* @return A tuple containing the x and y coordinate of the specified node.
*/
std::tuple<unsigned int, unsigned int> lbm_access::get_node_coordinates(unsigned int node_index)
{
    return std::make_tuple(node_index % HORIZONTAL_NODES, node_index / HORIZONTAL_NODES);
}

/**
 * @brief This function returns the distribution values of the node with the specified index using the specified access pattern.
 * 
 * @param source the distribution values will be read from this vector
 * @param node_index this is the index of the node in the domain
 * @param access this access function will be used
 * @return A vector containing the distribution values
 */
std::vector<double> lbm_access::get_distribution_values_of
(
    const std::vector<double> &source, 
    int node_index, 
    access_function access
)
{
    std::vector<double> dist_vals(9,0);
    for(auto direction = 0; direction < DIRECTION_COUNT; ++direction)
    {
        dist_vals[direction] = source[access(node_index, direction)];
    }
    return dist_vals;
}

/**
 * @brief This function sets all distribution values of the node with the specified index to the specified values 
 *        using the specified access pattern.
 * 
 * @param dist_vals a vector containing the values to which the distribution values shall be set
 * @param destination the distribution values will be written to this vector
 * @param node_index this is the index of the node in the domain
 * @param access this access function will be used
 */
void lbm_access::set_distribution_values_of
(
    const std::vector<double> &dist_vals, 
    std::vector<double> &destination, 
    int node_index, 
    access_function access
)
{
    for(auto direction = 0; direction < DIRECTION_COUNT; ++direction)
    {
        destination[access(node_index, direction)] = dist_vals[direction];
    }
}

/**
 * @brief Returns a vector containing all fluid segments in the specified domain. The order is as follows:
 *        All even numbers mark the beginning of a fluid node.
 *        All odd numbers state how many fluid tiles there are in a row, i.e. how many consecutive fluid tiles there are.
 *        For example, the sequence "3, 10" means that the nodes 3...12 are fluid nodes.
 * 
 * @param phase_space a vector containing the phase information for each node where "true" means solid
 * @return a vector containing the fluid segments in the explained arrangement
 */
std::vector<unsigned int> semi_direct_access::get_fluid_segments(const std::vector<bool> &node_phases)
{
    unsigned int index = 0;
    unsigned int consecution = 0;
    std::list<unsigned int> fluid_segments;

    while (index < size(node_phases))
    {
        if (!node_phases[index + consecution]) consecution++; // Hit fluid node
        else if (consecution > 0) // Hit solid node that marks the end of a line of consecutive fluid nodes
        {
            fluid_segments.push_back(index);
            fluid_segments.push_back(consecution);
            index += consecution;
            consecution = 0;
        }
        else index++; // Hit consecutive solid nodes
    }

    std::vector<unsigned int> result
        {
            std::make_move_iterator(begin(fluid_segments)), 
            std::make_move_iterator(end(fluid_segments))
        };

    return result;
}


