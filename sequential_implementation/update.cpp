#include "defines.hpp"
#include "boundaries.hpp"
#include "access.hpp"

/**
 * @brief This type represents a tuple containing all boundary nodes in the following order:
 *        Corner nodes (upper left, lower left, upper right, lower right) ->
 *        inlet nodes ->
 *        outlet nodes ->
 *        upper wall nodes ->
 *        lower wall nodes.
 *        In order to avoid unnecessary branching for comparably few nodes, these nodes will be treated separately by design.
 */
typedef std::tuple<std::array<int, 4>, 
                std::array<int, VERTICAL_NODES - 2>,
                std::array<int, VERTICAL_NODES - 2>,
                std::array<int, HORIZONTAL_NODES - 2>,
                std::array<int, HORIZONTAL_NODES - 2>
                > 
                node_tuple;

namespace update
{
    arr_of_dist_val collide_bgk(arr_of_dist_val f, velocity u, double density)
    {
        arr_of_dist_val result = maxwell_boltzmann_distribution(u, density);
        for(auto i = 0; i < DIRECTION_COUNT; ++i)
        {
            result[i] = -(1/RELAXATION_TIME) * (f[i] - result[i]);
        }
        return result;
    }

    void perform_two_step_stream(std::function<unsigned int(unsigned int, unsigned int)> access_function, std::array<double, TOTAL_NODE_COUNT>& nodes)
    {
        // Perform streaming step for all corner nodes
        /*
        nodes[access::get_node_index()]

        boundaries::upper_inlet_boundary_stream();
        */
    }

    /*
    void regular_outstream(std::function<unsigned int(unsigned int, unsigned int)> access_function, std::array<double, TOTAL_NODE_COUNT>& nodes)
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

            if(current_velocity_x == 1) // If facing rightward, increase column start offset by 1
            {
                column_start_offset = 1;
                column_end_offset = 0;
            }
            else if(current_velocity_x == 1) // If facing leftward, decrease column end offset by 1
            {
                column_start_offset = 0;
                column_end_offset = 1;
            }
            else // facing neigher leftward nor rightward
            {
                column_start_offset = 0;
                column_end_offset = 0;
            }
            
            
            
            for(auto row = row_start_offset; row < HORIZONTAL_NODES - 1 - row_end_offset; ++row)
            {
                for(auto column = column_start_offset; column < VERTICAL_NODES - 1 - column_end_offset; ++column)
                {

                 swap(
                            nodes, 
                            access::get_node_index(row, column), 
                            access::get_node_index(row + velocity_vectors[direction][0], column + velocity_vectors[direction][1])
                            );

                }
            }
        }
    }
    */
}
