#include "update.hpp"

arr_of_dist_val collision::collide_bgk(arr_of_dist_val f, velocity u, double density)
{
    arr_of_dist_val result = maxwell_boltzmann_distribution(u, density);
    for(auto i = 0; i < DIRECTION_COUNT; ++i)
    {
        result[i] = -(1/RELAXATION_TIME) * (f[i] - result[i]);
    }
    return result;
}

void stream::two_lattice::helper::two_lattice_corners(
    double &v_x, 
    double &v_y, 
    std::array<double, 15UL> &destination, 
    access_function &access_function, 
    std::array<double, 15UL> &source)
{
    for (auto corner : corner_specs)
    {
        for (auto direction : std::get<1>(corner))
        {
            v_x = velocity_vectors[direction][0];
            v_y = velocity_vectors[direction][1];
            destination[access_function(
                access::get_node_index(std::get<0>(std::get<0>(corner)) + v_x, std::get<1>(std::get<0>(corner)) + v_y),
                direction)] = source[access_function(access::get_node_index(std::get<0>(std::get<0>(corner)), std::get<1>(std::get<0>(corner))),
                                                        direction)];
        }
    }
}

void stream::two_lattice::helper::two_lattice_outlet(
    double &v_x, 
    double &v_y, 
    std::array<double, 15UL> &destination, 
    access_function &access_function, 
    std::array<double, 15UL> &source)
{
    for (auto row = 1; row < VERTICAL_NODES - 2; ++row)
    {
        for (auto direction : boundaries::neighbor_directions[boundaries::boundary_scenarios::outlet])
        {
            v_x = velocity_vectors[direction][0];
            v_y = velocity_vectors[direction][1];
            destination[access_function(access::get_node_index(HORIZONTAL_NODES - 1 + v_x, row + v_y), direction)] = source[access_function(access::get_node_index(HORIZONTAL_NODES - 1, row), direction)];
        }
    }
}

void stream::two_lattice::helper::two_lattice_inlet(
    double &v_x, 
    double &v_y, 
    std::array<double, 15UL> &destination, 
    access_function &access_function, 
    std::array<double, 15UL> &source)
{
    for (auto row = 1; row < VERTICAL_NODES - 2; ++row)
    {
        for (auto direction : boundaries::neighbor_directions[boundaries::boundary_scenarios::inlet])
        {
            v_x = velocity_vectors[direction][0];
            v_y = velocity_vectors[direction][1];
            destination[access_function(access::get_node_index(v_x, row + v_y), direction)] = source[access_function(access::get_node_index(0, row), direction)];
        }
    }
}

void stream::two_lattice::helper::two_lattice_walldown(
    double &v_x, 
    double &v_y, 
    std::array<double, 15UL> &destination, 
    access_function &access_function, 
    std::array<double, 15UL> &source)
{
    for (auto column = 1; column < HORIZONTAL_NODES; ++column)
    {
        for (auto direction : boundaries::neighbor_directions[boundaries::boundary_scenarios::wall_down])
        {
            v_x = velocity_vectors[direction][0];
            v_y = velocity_vectors[direction][1];
            destination[access_function(access::get_node_index(column + v_x, VERTICAL_NODES - 1 + v_y), direction)] = source[access_function(access::get_node_index(column, VERTICAL_NODES - 1), direction)];
        }
    }
}

void stream::two_lattice::helper::two_lattice_wallup(
    double &v_x, 
    double &v_y,
    std::array<double, 15UL> &destination, 
    access_function &access_function, 
    std::array<double, 15UL> &source)
{
    for (auto column = 1; column < HORIZONTAL_NODES; ++column)
    {
        for (auto direction : boundaries::neighbor_directions[boundaries::boundary_scenarios::wall_up])
        {
            v_x = velocity_vectors[direction][0];
            v_y = velocity_vectors[direction][1];
            destination[access_function(access::get_node_index(column + v_x, v_y), direction)] = source[access_function(access::get_node_index(column, 0), direction)];
        }
    }
}

void stream::two_lattice::helper::two_lattice_regular(
    double &v_x, 
    double &v_y, 
    std::array<double, 15UL> &destination, 
    access_function &access_function, 
    std::array<double, 15UL> &source)
{
    for (auto row = 1; row < VERTICAL_NODES - 2; ++row)
    {
        for (auto column = 1; column < HORIZONTAL_NODES - 2; ++column)
        {
            for (auto direction : general_stream_directions)
            {
                v_x = velocity_vectors[direction][0];
                v_y = velocity_vectors[direction][1];
                destination[access_function(access::get_node_index(column + v_x, row + v_y), direction)] = source[access_function(access::get_node_index(column, row), direction)];
            }
        }
    }
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

