#include "defines.hpp"
#include "boundaries.hpp"
#include "access.hpp"
#include <array>

/**
 * @brief This type represents a tuple containing all boundary nodes in the following order:
 *        Corner nodes (upper left, lower left, upper right, lower right) ->
 *        inlet nodes ->
 *        outlet nodes ->
 *        upper wall nodes ->
 *        lower wall nodes.
 *        In order to avoid unnecessary branching for comparably few nodes, these nodes will be treated separately by design.
 */
typedef std::tuple<
                std::array<int, 4>, 
                std::array<int, VERTICAL_NODES - 2>,
                std::array<int, VERTICAL_NODES - 2>,
                std::array<int, HORIZONTAL_NODES - 2>,
                std::array<int, HORIZONTAL_NODES - 2>
                > 
                node_tuple;

/**
 * @brief This type stands for an access function. Node values can be stored in different layout and via this function,
 *        the corresponding access scheme can be specified.
 */
typedef std::function<unsigned int(unsigned int, unsigned int)> access_function;

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
     * @return arr_of_dist_val 
     */
    arr_of_dist_val collide_bgk(arr_of_dist_val f, velocity u, double density);
}

/**
 * @brief This namespace contains a namespace for each stream algorithm. 
 *        Namespaces named after an algorithm provide functions that are meant for access from outside.
 *        Any namespaces declared as "helper" are internal only and not intended for outside access.
 * 
 */
namespace stream
{
    /**
     * @brief Contains the directions that are considered during streaming for regular (non-border) nodes.
     * 
     */
    std::array<unsigned int, DIRECTION_COUNT - 1> general_stream_directions{0,1,2,3,5,6,7,8};

    /**
     * @brief This tuple contains the following information that is necessary for the streaming of border node values
     * 
     */
    std::array<std::tuple<std::tuple<unsigned int, unsigned int>, std::list<int>>, 4> corner_specs{
        std::make_tuple(std::make_tuple(0,0), boundaries::neighbor_directions[boundaries::boundary_scenarios::lower_inlet]),
        std::make_tuple(std::make_tuple(0,VERTICAL_NODES - 1), boundaries::neighbor_directions[boundaries::boundary_scenarios::upper_inlet]),
        std::make_tuple(std::make_tuple(HORIZONTAL_NODES - 1,0), boundaries::neighbor_directions[boundaries::boundary_scenarios::lower_outlet]),
        std::make_tuple(std::make_tuple(HORIZONTAL_NODES - 1,VERTICAL_NODES - 1), boundaries::neighbor_directions[boundaries::boundary_scenarios::upper_outlet])
    };

    /**
     * @brief This namespace provides a method to perform a streaming step using the two-lattice algorithm.
     *        It also includes various helper functions in a separate namespace which are not meant for outside access.
     * 
     */
    namespace two_lattice
    {
        namespace helper
        {
            // Perform streaming for all regular nodes
            void two_lattice_regular(
                double &v_x, double &v_y, std::array<double, 15UL> &destination, 
                access_function &access_function, 
                std::array<double, 15UL> &source
                );

            // Perform streaming step for all upper wall nodes   
            void two_lattice_wallup(
                double &v_x, double &v_y, 
                std::array<double, 15UL> &destination, 
                access_function &access_function, 
                std::array<double, 15UL> &source);

            // Perform streaming step for all lower-wall nodes   
            void two_lattice_walldown(
                double &v_x, double &v_y, 
                std::array<double, 15UL> &destination, 
                access_function &access_function, 
                std::array<double, 15UL> &source);

            // Perform streaming step for all inlet-only nodes    
            void two_lattice_inlet(
                double &v_x, double &v_y, 
                std::array<double, 15UL> &destination, 
                access_function &access_function, 
                std::array<double, 15UL> &source);

            // Perform streaming step for all outlet-only nodes    
            void two_lattice_outlet(
                double &v_x, double &v_y, 
                std::array<double, 15UL> &destination, 
                access_function &access_function, 
                std::array<double, 15UL> &source);

            // Perform streaming step for all corner nodes    
            void two_lattice_corners(
                double &v_x, double &v_y, 
                std::array<double, 15UL> &destination, 
                access_function &access_function, 
                std::array<double, 15UL> &source);
        }


        /**
         * @brief 
         * 
         * @param access_function 
         * @param nodes_0 
         * @param nodes_1 
         * @param timestep 
         */
        void perform_two_lattice_stream(
            access_function access_function,
            std::array<double, TOTAL_NODE_COUNT> &nodes_0,
            std::array<double, TOTAL_NODE_COUNT> &nodes_1,
            unsigned int timestep)
        {
            double v_x = 0;
            double v_y = 0;
            std::array<double, TOTAL_NODE_COUNT>& source = (timestep % 2) ? nodes_0 : nodes_1;
            std::array<double, TOTAL_NODE_COUNT>& destination = (timestep % 2) ? nodes_1 : nodes_0;

            // Perform streaming
            helper::two_lattice_regular(v_x, v_y, destination, access_function, source);
            helper::two_lattice_wallup(v_x, v_y, destination, access_function, source);
            helper::two_lattice_walldown(v_x, v_y, destination, access_function, source);
            helper::two_lattice_inlet(v_x, v_y, destination, access_function, source);
            helper::two_lattice_outlet(v_x, v_y, destination, access_function, source);
            helper::two_lattice_corners(v_x, v_y, destination, access_function, source);

            for (int direction; direction < DIRECTION_COUNT; ++direction)
            {
                // Upper inlet
                /*
                destination[
                    access_function(
                        access::get_node_index(std::get<0>(std::get<0>(corner)), std::get<1>(std::get<0>(corner))),
                        direction)
                ] = boundaries::
                */

                // Lower inlet

                // Upper outlet

                // Lower outlet

                // Upper wall

                // Lower wall

                // inlet

                // outlet

                for (auto direction : boundaries::neighbor_directions[boundaries::boundary_scenarios::outlet])
                {
                    for (auto row = 1; row < VERTICAL_NODES - 2; ++row)
                    {
                        v_x = velocity_vectors[direction][0];
                        v_y = velocity_vectors[direction][1];
                        destination[access_function(access::get_node_index(HORIZONTAL_NODES - 1, row), direction)] = boundaries::outlet_boundary_stream()
                    }

    


                }

            }
        }
    }

}