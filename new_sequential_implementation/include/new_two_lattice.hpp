#ifndef TWO_LATTICE_SEQUENTIAL_HPP
#define TWO_LATTICE_SEQUENTIAL_HPP

#include <vector>
#include <set>
#include "defines.hpp"
#include "access.hpp"
#include <iostream>
#include "new_collision.hpp"
#include "utils.hpp"
#include "boundaries.hpp"

namespace two_lattice_sequential
{

    /**
     * @brief Performs the combined streaming and collision step for all fluid nodes within the simulation domain.
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
        std::vector<unsigned int> &fluid_nodes,
        border_swap_information &bsi,
        std::vector<double> &source, 
        std::vector<double> &destination,    
        access_function access_function
    );

    /**
     * @brief Performs the sequential two-lattice algorithm for the specified number of iterations.
     * 
     * @param fluid_nodes A vector containing the indices of all fluid nodes in the domain
     * @param boundary_nodes A vector containing the indices of all fluid boundary nodes in the domain
     * @param values_0 source for even time steps and destination for odd time steps
     * @param values_1 source for odd time steps and destination for even time steps
     * @param access_function the access function according to which distribution values are to be accessed
     * @param iterations this many iterations will be performed
     * @param data the simulation data tuples will be placed in this vector (assumed to be pre-initialized)
     */
    void run
    (  
        std::vector<unsigned int> &fluid_nodes,       
        border_swap_information &boundary_nodes,
        std::vector<double> &values_0, 
        std::vector<double> &values_1,   
        access_function access_function,
        unsigned int iterations,
        std::vector<sim_data_tuple> &data
    );

    /**
     * @brief Determines the remaining streaming option for a node based on the specified border 
     *        information vector.
     * 
     * @param current_border_info an entry of a border_swap_information object
     * @return a set containing all remaining streaming directions
     */
    inline std::set<unsigned int> determine_streaming_directions
    (
        std::vector<unsigned int> &current_border_info
    )
    {
        std::set<unsigned int> remaining_dirs = {streaming_directions.begin(), streaming_directions.end()};
        std::set<unsigned int> bounce_back_dirs = bounce_back::determine_bounce_back_directions(current_border_info);
        // if(current_border_info[0] == 19)
        // {
        //     std::cout << "Bounce back dirs for node 19 are ";
        //     to_console::print_set(bounce_back_dirs); 
        // }
        
        for (auto i : bounce_back_dirs)
        {
            // if(current_border_info[0] == 19)
            // {
            //     std::cout << "Removing direction " << i << std::endl;
            // }
            remaining_dirs.erase(i);
        }
        unsigned int x = std::get<0>(access::get_node_coordinates(current_border_info[0]));
        unsigned int y = std::get<1>(access::get_node_coordinates(current_border_info[0]));
        // if(current_border_info[0] == 19)
        // {
        //     std::cout << "x coordinate is " << std::get<0>(access::get_node_coordinates(current_border_info[0])) << std::endl;
        // }
        if(x == 1)
        {
            if(y == 1) remaining_dirs.insert({2,5});
            else if(y == (VERTICAL_NODES - 2)) remaining_dirs.insert({5,8});
            else remaining_dirs.insert({2,5,8});
            
        }
        else if(x ==(HORIZONTAL_NODES - 2))
        {
            if(y == 1) remaining_dirs.insert({0,3});
            else if(y == (VERTICAL_NODES - 2)) remaining_dirs.insert({3,6});
            else remaining_dirs.insert({0,3,6});
        }
        if(current_border_info[0] == 19)
        {
            std::cout << "Streaming dirs for node 19 are ";
            to_console::print_set(remaining_dirs); 
        }
        return remaining_dirs;  
        // std::set<unsigned int> remaining_dirs = {streaming_directions.begin(), streaming_directions.end()};
        // std::set<unsigned int> bounce_back_dirs = bounce_back::determine_bounce_back_directions(current_border_info);
        // //std::cout << "Received current border info ";
        // //to_console::print_vector(current_border_info, 10);
        // //std::cout << "Bounce back dirs are ";
        // //to_console::print_set(bounce_back_dirs);
        // for (auto i : bounce_back_dirs)
        // {
        //     remaining_dirs.erase(invert_direction(i));
        // }
        // //std::cout << "Returning remaining dirs ";
        // //to_console::print_set(remaining_dirs);
        // return remaining_dirs;  
    }

    /**
     * @brief Performs the steaming step in the specified directions for the fluid node with 
     *        the specified index.
     * 
     * @param source distribution values will be taken from this vector
     * @param destination distribution values will be rearranged in this vector
     * @param access_function function that will be used to access the distribution values
     * @param fluid_node the index of the node for which the streaming step is performed
     * @param directions a set specifying in which directions streaming will be executed
     */
    inline void tl_stream
    (
        std::vector<double> &source,
        std::vector<double> &destination, 
        access_function &access_function, 
        unsigned int fluid_node, 
        std::set<unsigned int> &directions
    )
    {
        std::cout << "Executing streaming step for node " << fluid_node << std::endl;
        std::cout << "Got directions ";
        to_console::print_set(directions);
        for(auto direction : directions) 
        {
            std::cout << "\t performing stream (node = " << fluid_node << ", dir = " << direction << ") := " << "(node = " << access::get_neighbor(fluid_node, invert_direction(direction)) << ", dir = " << direction << ")" << std::endl; 
            destination[access_function(fluid_node, direction)] =
            source[
                access_function(
                    access::get_neighbor(fluid_node, invert_direction(direction)), 
                    direction)];
        }
    }

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
        std::vector<double> &source,
        std::vector<double> &destination, 
        access_function &access_function, 
        unsigned int fluid_node
    )
    {
        for (auto direction : streaming_directions)
        {
            destination[access_function(fluid_node, direction)] =
                source[
                    access_function(
                        access::get_neighbor(fluid_node, invert_direction(direction)), 
                        direction)];
        }

    }

    /**
     * @brief Performs the collision step for the fluid node with the specified index.
     * 
     * @param destination the updated distribution values will be written to this vector
     * @param fluid_node the index of the fluid node
     * @param access_function the function used to access the distribution values
     * @param velocities a vector containing the velocities of ALL nodes within the lattice
     * @param densities a vector containing the densities of ALL nodes within the lattice
     */
    inline void tl_collision
    (
        std::vector<double> &destination, 
        unsigned int fluid_node, 
        std::vector<double> &distribution_values,
        access_function &access_function, 
        std::vector<velocity> &velocities, 
        std::vector<double> &densities
    )
    {
        //std::cout << "!!!!!!!!!!!!!!!!!!!!!!! Accessin collide_bgk with vals of length " << distribution_values.size() << std::endl;
        std::vector<double >vals = collision::collide_bgk(distribution_values, velocities[fluid_node], densities[fluid_node]);
        access::set_all_distribution_values(vals, destination, fluid_node, access_function);
    }
}

#endif