#ifndef NEW_COLLISION_HPP
#define NEW_COLLISION_HPP
#include "defines.hpp"

/**
 * @brief This namespace contains all functions necessary for the collision step.
 */
namespace collision
{
    /**
     * @brief Performs the collision step for a node with the specified distribution values, velocity and density.
     * 
     * @param values a vector containing the distribution values of this node.
     * @param u the flow velocity at this node
     * @param density the density at this node
     * @return a vector containing the updated distribution values.
     */
    std::vector<double> collide_bgk
    (
        const std::vector<double> &values, 
        velocity &u, 
        double density
    );
    
    /**
     * @brief Performs the collision step for all fluid nodes.
     *        This function is intended to be used in the two-step algorithm as the streaming and collision steps cannot be fused there.
     * 
     * @param fluid_nodes a vector containing the node indices of all fluid nodes
     * @param values a vector containing the distribution values of all nodes
     * @param all_velocities a vector containing the velocity values of all fluid nodes
     * @param all_densities a vector containing the density values of all fluid nodes
     * @param access this function is used to access the distribution values
     */
    void collide_all_bgk
    (
        const std::vector<unsigned int> &fluid_nodes,
        std::vector<double> &values, 
        std::vector<velocity> &all_velocities, 
        std::vector<double> &all_densities,
        access_function access
    );
}
#endif