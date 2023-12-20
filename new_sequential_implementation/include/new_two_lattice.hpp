#ifndef TWO_LATTICE_SEQUENTIAL_HPP
#define TWO_LATTICE_SEQUENTIAL_HPP

#include <vector>
#include <set>
#include "defines.hpp"
#include "access.hpp"

namespace two_lattice_sequential
{

    /**
     * @brief Performs the combined streaming-and-collision step for all nodes within the simulation domain.
     * 
     * @param time_step current iteration of the simulation
     * @param fluid_nodes A vector containing the indices of all fluid nodes in the domain
     * @param boundary_nodes A vector containing the indices of all fluid boundary nodes in the domain
     * @param values_0 source for even time steps and destination for odd time steps
     * @param values_1 source for odd time steps and destination for even time steps
     * @param access_function the access function according to which the values are to be accessed
     * @return a tuple containing vectors of all flow velocities and density values for a fixed time step.
     */
    void perform_tl_stream_and_collide
    (
        std::vector<unsigned int> &fluid_nodes,
        border_swap_information &bsi,
        std::vector<double> &source, 
        std::vector<double> &destination,    
        access_function access_function,
        sim_data_vector &sim_data
    );

    /**
     * @brief Performs the sequential two-lattice algorithm for the specified number of iterations.
     * 
     * @param fluid_nodes A vector containing the indices of all fluid nodes in the domain
     * @param boundary_nodes A vector containing the indices of all fluid boundary nodes in the domain
     * @param values_0 source for even time steps and destination for odd time steps
     * @param values_1 source for odd time steps and destination for even time steps
     * @param access_function the access function according to which the values are to be accessed
     * @param iterations this many iterations will be performed
     * @return a vector of tuples containing all flow velocities and density values for all time steps
     */
    void run
    (  
        std::vector<unsigned int> &fluid_nodes,       
        border_swap_information &boundary_nodes,
        std::vector<double> &values_0, 
        std::vector<double> &values_1,   
        access_function access_function,
        unsigned int iterations,
        sim_data_vector &data
    );

    /**
     * @brief 
     * 
     * @param current_border_info 
     * @return std::set<unsigned int> 
     */
    inline std::set<unsigned int> determine_remaining_directions
    (
        std::vector<unsigned int> &current_border_info
    )
    {
        std::set<unsigned int> remaining_dirs = {streaming_directions.begin(), streaming_directions.end()};
        for (auto i = current_border_info.begin() + 1; i < current_border_info.end(); ++i)
        {
            remaining_dirs.erase(*i);
        }
        return remaining_dirs;  
    }

    /**
     * @brief 
     * 
     * @param source 
     * @param destination 
     * @param access_function 
     * @param fluid_node 
     * @param direction 
     */
    inline void tl_stream
    (
        std::vector<double> &source,
        std::vector<double> &destination, 
        access_function &access_function, 
        unsigned int fluid_node, 
        int direction 
    )
    {
        destination[access_function(fluid_node, direction)] =
            source[access_function(access::get_neighbor(fluid_node, direction), invert_direction(direction))];
    }

    /**
     * @brief 
     * 
     * @param destination 
     * @param fluid_node 
     * @param access_function 
     * @param velocities 
     * @param densities 
     */
    void tl_collision
    (
        std::vector<double> &destination, 
        unsigned int fluid_node, 
        access_function &access_function, 
        std::vector<velocity> &velocities, 
        std::vector<double> &densities
    );

}

#endif