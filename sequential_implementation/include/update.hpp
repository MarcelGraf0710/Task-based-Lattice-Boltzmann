#ifndef UPDATE_HPP
#define UPDATE_HPP
#include "defines.hpp"
#include "boundaries.hpp"
#include "access.hpp"
#include <array>

/**
 * @brief This namespace contains all functions necessary for the collision step.
 */
namespace collision
{
    /**
     * @brief Performs the collision step using the specified array of distribution values and the parameters necessary for determining
     *        the equilibrium distribution.
     * 
     * @param f an array containing the current distribution values
     * @param u the velocity at this node
     * @param density the density at this node
     * @return an array containing all distribution values after the collision step
     */
    vec_of_dist_val collide_bgk(vec_of_dist_val f, velocity u, double density);

    /**
     * @brief 
     * 
     * @param all_distribution_values 
     * @param all_velocities 
     * @param access 
     * @param density 
     */
    void perform_collision_step(std::vector<double> all_distribution_values, std::vector<velocity> all_velocities, access_function access, double density = INFINITY);
}

/**
 * @brief This namespace contains a namespace for each stream algorithm. 
 *        Namespaces named after an algorithm provide functions that are meant for access from outside.
 *        Any namespaces declared as "helper" are internal only and not intended for outside access.
 */
namespace stream
{
    /**
     * @brief Contains the directions that are considered during streaming for regular (non-border) nodes.
     */
    std::array<unsigned int, DIRECTION_COUNT - 1> general_stream_directions{0,1,2,3,5,6,7,8};

    /**
     * @brief This tuple contains the following information that is necessary for the streaming of border node values
     */
    std::array<std::tuple<std::tuple<unsigned int, unsigned int>, std::list<int>>, 4> corner_specs{
        std::make_tuple(std::make_tuple(0,0), std::list{5,7,8}),
        std::make_tuple(std::make_tuple(0,VERTICAL_NODES - 1), std::list{1,2,5}),
        std::make_tuple(std::make_tuple(HORIZONTAL_NODES - 1,0), std::list{3,6,7}),
        std::make_tuple(std::make_tuple(HORIZONTAL_NODES - 1,VERTICAL_NODES - 1), std::list{0,1,3})
    };

    /**
     * @brief This namespace provides a method to perform a streaming step using the two-lattice algorithm.
     *        It also includes various helper functions in a separate namespace which are not meant for outside access.
     */
    namespace two_lattice
    {
        namespace helper
        {
            // Perform streaming for all regular nodes
            void two_lattice_regular(
                all_distributions &source,
                all_distributions &destination, 
                access_function &access_function
                );

            // Perform streaming step for all upper wall nodes   
            void two_lattice_wallup(
                all_distributions &source,
                all_distributions &destination, 
                access_function &access_function
                );

            // Perform streaming step for all lower-wall nodes   
            void two_lattice_walldown(
                all_distributions &source,
                all_distributions &destination, 
                access_function &access_function
                );

            // Perform streaming step for all inlet-only nodes    
            void two_lattice_inlet(
                all_distributions &source,
                all_distributions &destination, 
                access_function &access_function
                );

            // Perform streaming step for all outlet-only nodes    
            void two_lattice_outlet(
                all_distributions &source,
                all_distributions &destination, 
                access_function &access_function
                );

            // Perform streaming step for all corner nodes    
            void two_lattice_corners(
                all_distributions &source,
                all_distributions &destination, 
                access_function &access_function
                );

            // Updates the distribution values that are missing after streaming
            void perform_boundary_update(
                all_distributions &source, 
                all_distributions &destination, 
                access_function &access_function);
        }

        /**
         * @brief Performs one streaming step for all nodes using the two-lattice algorithm.
         *        During the process, all necessary post-streaming updates of border nodes are performed such that
         *        the collision step can be executed next without any further redo.
         * 
         * @param access_function this access function determines the pattern in which distribution values are accessed within the arrays containing 
         *                        the distribution values
         * @param nodes_0 this node stores all node distribution values and acts as a source during even time steps
         * @param nodes_1 this node stores all node distribution values and acts as a source during odd time steps
         * @param timestep the module of this timestep is used to determine the source and destination roles of the two specified arrays
         */
        void perform_two_lattice_stream(
            access_function access_function,
            all_distributions &source,
            all_distributions &destination
            )
        {
            helper::two_lattice_regular(source, destination, access_function);
            helper::two_lattice_wallup(source, destination, access_function);
            helper::two_lattice_walldown(source, destination, access_function);
            helper::two_lattice_inlet(source, destination, access_function);
            helper::two_lattice_outlet(source, destination, access_function);
            helper::two_lattice_corners(source, destination, access_function);
            helper::perform_boundary_update(source, destination, access_function);
        }
        
    }

    /**
     * @brief This namespace provides a method to perform a streaming step using the two-step algorithm.
     *        It also includes various helper functions in a separate namespace which are not meant for outside access.
     */
    namespace two_step
    {
        /**
         * @brief 
         * 
         * @param access_function 
         * @param source 
         */
        void perform_two_lattice_stream(
            access_function access_function,
            all_distributions &source
        )
        {
            unsigned int row_start_offset = 0;
            unsigned int row_end_offset = 0;
            unsigned int column_start_offset = 0;
            unsigned int column_end_offset = 0;
            double current_velocity_x = 0;
            double current_velocity_y = 0;
            for(auto direction = 0; direction < DIRECTION_COUNT; ++direction)
            {
                current_velocity_x = velocity_vectors[direction][0];
                current_velocity_y = velocity_vectors[direction][1];
                if(current_velocity_y == -1) // If facing downward, increase row start offset by 1
                {
                    row_start_offset = 1;
                    row_end_offset = 0;
                }
                else if(current_velocity_y == 1) // If facing upward, decrease row end offset by 1 
                {
                    row_start_offset = 0;
                    row_end_offset = 1;
                }
                 else // facing neigher upward nor downward
                {
                    row_start_offset = 0;
                    row_end_offset = 0;
                }

                if(current_velocity_x == -1) // If facing rightward, increase column start offset by 1
                {
                    column_start_offset = 1;
                    column_end_offset = 0;
                }
                else if(current_velocity_x == 1) // If facing leftward, decrease column end offset by 1 
                {
                    column_start_offset = 0;
                    column_end_offset = 1;
                }
                 else // facing neigher rightward nor leftward
                {
                    column_start_offset = 0;
                    column_end_offset = 0;
                }
            }
        }


    }
}

#endif
