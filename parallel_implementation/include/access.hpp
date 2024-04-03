#ifndef ACCESS_HPP
#define ACCESS_HPP
#include "defines.hpp"

/**
 * @brief This namespace contains functions that map input values to array index accesses.
 */
namespace lbm_access
{
   /**
    * @brief Retrieves the coordinates of the node with the specified node index.
    * @return A tuple containing the x and y coordinate of the specified node.
    */
    std::tuple<unsigned int, unsigned int> get_node_coordinates(unsigned int node_index);

    /**
     * @brief Returns the index of the neighbor that is reached when moving in the specified direction.
     * 
     * @param node_index the index of the current node
     * @param direction the direction of movement
     * @return the node index of the neighbor
     */
    inline unsigned int get_neighbor(unsigned int node_index, unsigned int direction)
    {
        int y_offset = direction / 3 - 1; // -1 for {0,1,2}, 0 for {3,4,5}, 1 for {6,7,8}
        int x_offset = direction - (3 * y_offset + 4); //-1 for {0,3,6}, 0 for {1,4,7}, 1 for {2,5,8}
        return node_index + y_offset * HORIZONTAL_NODES + x_offset;
    }

    /**
     * @brief Returns the array index of the collision layout
     * 
     * @param node the node in the simulation domain
     * @param direction the direction of the velocity vector
     * @return the index of the array storing the distribution values 
     */
    inline unsigned int collision(unsigned int node, unsigned int direction)
    {
        return DIRECTION_COUNT * node + direction;        
    }

    /**
     * @brief Returns the array index of the stream layout
     * 
     * @param node the node in the simulation domain
     * @param direction the direction of the velocity vector
     * @return the index of the array storing the distribution values  
     */
    inline unsigned int stream(unsigned int node, unsigned int direction)
    {
        return TOTAL_NODE_COUNT * direction + node;
    }

    /**
     * @brief Returns the array index of the bundle layout
     * 
     * @param node the node in the simulation domain
     * @param direction the direction of the velocity vector
     * @return the index of the array storing the distribution values  
     */
    inline unsigned int bundle(unsigned int node, unsigned int direction)
    {
        return 3 * (direction / 3) * TOTAL_NODE_COUNT + (direction % 3) + 3 * node; 
    }

    /**
     * @brief Returns the index the desired node has within the array that stores it. 
     *        The origin lies at the lower left corner and enumeration is row-major.
     * 
     * @param x x coordinate
     * @param y y coordinate
     * @return the index of the desired note
     */
    inline unsigned int get_node_index(unsigned int x, unsigned int y)
    {
        return x + y * HORIZONTAL_NODES;
    }

    /**
     * @brief This function returns the distribution values of the node with the specified index using the specified access pattern.
     * 
     * @param source the distribution values will be read from this vector
     * @param node_index this is the index of the node in the domain
     * @param access this access function will be used
     * @return A vector containing the distribution values
     */
    std::vector<double> get_distribution_values_of
    (
        const std::vector<double> &source, 
        int node_index, 
        access_function access
    );

    /**
     * @brief This function sets the distribution values of the node with the specified index to the specified values 
     *        using the specified access pattern.
     * 
     * @param dist_vals a vector containing the values to which the distribution values shall be set
     * @param destination the distribution values will be written to this vector
     * @param node_index this is the index of the node in the domain
     * @param access this access function will be used
     */
    void set_distribution_values_of
    (
        const std::vector<double> &dist_vals, 
        std::vector<double> &destination, 
        int node_index, 
        access_function access
    );
}

/**
 * @brief This namespace contains utility functions for semi-direct access schemes.
 *        Currently, only direct access is used. It is possible, however, to implement semi-direct access.
 */
namespace semi_direct_access
{
    /**
     * @brief Returns a vector containing all fluid segments in the specified domain. The order is as follows:
     *        All even numbers mark the beginning of a fluid node.
     *        All odd numbers state how many fluid tiles there are in a row, i.e. how many consecutive fluid tiles there are.
     *        For example, the sequence "3, 10" means that the nodes 3...12 are fluid nodes.
     * 
     * @param phase_space a vector containing the phase information for each node where "true" means solid
     * @return a vector containing the fluid segments in the explained arrangement
     */
    std::vector<unsigned int> get_fluid_segments(const std::vector<bool> &node_phases);
}

namespace buffered_access 
{

}

#endif