#include "../include/defines.hpp" 
#ifndef ACCESS_HPP
#define ACCESS_HPP

/**
 * @brief This namespace contains functions that map input values to array index accesses.
 */
namespace access
{
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
        return DIRECTION_COUNT * direction + node;
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
        return (direction / 3) * TOTAL_NODE_COUNT + (direction % 3); 
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
     * @brief This function gets all distribution values of the node with the specified index using the specified access pattern.
     * 
     * @param source the distribution values will be read from this vector
     * @param node_index this is the index of the node in the domain
     * @param access this access function will be used
     * @return All distribution values
     */
    std::vector<double> get_all_distribution_values(std::vector<double> &source, int node_index, access_function access);

    /**
     * @brief This function sets all distribution values of the node with the specified index to the specified values using the specified access pattern.
     * 
     * @param dist_vals a vector containing the values to which the distribution values shall be set
     * @param destination the distribution values will be written to this vector
     * @param node_index this is the index of the node in the domain
     * @param access this access function will be used
     */
    void set_all_distribution_values(std::vector<double> &dist_vals, std::vector<double> &destination, int node_index, access_function access);
}

/**
 * @brief This namespace contains utility functions for semi-direct access schemes.
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
    std::vector<unsigned int> get_fluid_segments(std::vector<bool> node_phases);
}

#endif