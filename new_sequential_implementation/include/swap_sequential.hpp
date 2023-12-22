#ifndef SWAP_SEQUENTIAL_HPP
#define SWAP_SEQUENTIAL_HPP

#include <vector>
#include "defines.hpp"

extern std::vector<unsigned int> general_swap_partners{5, 6, 7, 8};

namespace swap_sequential
{
    /**
     * @brief Sets up all distribution function values for boundary nodes
     * 
     * @param values the vector containing the distribution values of all nodes
     * @param access_function the access function according to which the values are to be accessed
     */
    void swap_boundary_initialization
    (
        std::vector<double> values, 
        access_function access_function
    );

    /**
     * @brief Utility function to determine all directions of non-solid swap partner nodes in the case of direct addressing.
     * 
     * @param current_node the node that is currently processed
     * @param all_nodes a vector containing all nodes in the simulation domain
     * @param node_phases a vector containing the phase information of all nodes in the simulation domain
     * @return a vector containing all directions pointing to non-solid swap partner nodes
     */
    std::vector<unsigned int> get_swap_partner_indices
    (
        unsigned int current_node,
        std::vector<int> all_nodes,
        std::vector<bool> node_phases 
    );

    /**
     * @brief Utility function to determine all directions of non-solid swap partner nodes in the case of 
     *        semi-direct addressing as proposed by Mattila et al.
     * 
     * @param current_node the node that is currently processed
     * @param phase_vector a vector containing the phase information of all nodes in the simulation domain
     * @return a vector containing all directions pointing to non-solid swap partner nodes
     */
    std::vector<unsigned int> get_swap_partner_indices
    (
        unsigned int current_node,
        std::vector<int> phase_vector
    );

    /**
     * @brief Swaps the specified elements of the specified vector of distribution functions.
     * 
     * @param values the vector containing the distribution values of all nodes
     * @param access_function the access function according to which the values are to be accessed
     * @param index_0 first swap index
     * @param index_1 second swap index
     */
    void swap_values
    (
        std::vector<double> values, 
        access_function access_function,
        unsigned int index_0,
        unsigned int index_1
    );

    /**
     * @brief Corrects the order of the distribution functions before collision.
     * 
     * @param values the vector containing the distribution values of all nodes
     * @param access_function the access function according to which the values are to be accessed
     */
    void correct_distribution_order
    (        
        std::vector<double> values, 
        access_function access_function
    );

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
        std::vector<unsigned int> fluid_nodes,       
        std::vector<double> values, 
        access_function access_function,
        unsigned int iterations
    );
}

#endif