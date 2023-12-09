#include "defines.hpp" 

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
     * @brief Get the all distribution values object
     * 
     * @param source 
     * @param node_index 
     * @param access 
     * @return std::vector<double> 
     */
    std::vector<double> get_all_distribution_values(std::vector<double> source, int node_index, access_function access)
    {
        std::vector<double> dist_vals;
        dist_vals.reserve(9);
        for(auto direction = 0; direction < DIRECTION_COUNT; ++direction)
        {
            dist_vals[direction] = source[access(node_index, direction)];
        }
        return dist_vals;
    }

    /**
     * @brief Set the all distribution values object
     * 
     * @param destination 
     * @param node_index 
     * @param access 
     * @return std::vector<double> 
     */
    void set_all_distribution_values(std::vector<double> dist_vals, std::vector<double> destination, int node_index, access_function access)
    {
        for(auto direction = 0; direction < DIRECTION_COUNT; ++direction)
        {
            destination[access(node_index, direction)] = dist_vals[direction];
        }
    }
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
    std::vector<unsigned int> get_fluid_segments(std::vector<bool> node_phases)
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
}