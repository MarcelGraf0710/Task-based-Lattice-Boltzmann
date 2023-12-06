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
    all_distributions &source,
    all_distributions &destination, 
    access_function &access_function
    )
{
    for (auto corner : corner_specs)
    {
        for (auto direction : std::get<1>(corner))
        {
            destination[access_function(
                access::get_node_index(std::get<0>(std::get<0>(corner)) + velocity_vectors[direction][0], std::get<1>(std::get<0>(corner)) + velocity_vectors[direction][1]),
                direction)] = source[access_function(access::get_node_index(std::get<0>(std::get<0>(corner)), std::get<1>(std::get<0>(corner))),
                                                        direction)];
        }
    }
}

void stream::two_lattice::helper::two_lattice_outlet(
    all_distributions &source,
    all_distributions &destination, 
    access_function &access_function
    )
{
    for (auto row = 1; row < VERTICAL_NODES - 2; ++row)
    {
        for (auto direction : boundaries::neighbor_directions[boundaries::boundary_scenarios::outlet])
        {
            destination[access_function(access::get_node_index(HORIZONTAL_NODES - 1 + velocity_vectors[direction][0], row + velocity_vectors[direction][1]), direction)] = source[access_function(access::get_node_index(HORIZONTAL_NODES - 1, row), direction)];
        }
    }
}

void stream::two_lattice::helper::two_lattice_inlet(
    all_distributions &source,
    all_distributions &destination, 
    access_function &access_function
    )
{
    for (auto row = 1; row < VERTICAL_NODES - 2; ++row)
    {
        for (auto direction : boundaries::neighbor_directions[boundaries::boundary_scenarios::inlet])
        {
            destination[access_function(access::get_node_index(velocity_vectors[direction][0], row + velocity_vectors[direction][1]), direction)] = source[access_function(access::get_node_index(0, row), direction)];
        }
    }
}

void stream::two_lattice::helper::two_lattice_walldown(
    all_distributions &source,
    all_distributions &destination, 
    access_function &access_function
    )
{
    for (auto column = 1; column < HORIZONTAL_NODES; ++column)
    {
        for (auto direction : boundaries::neighbor_directions[boundaries::boundary_scenarios::wall_down])
        {
            destination[access_function(access::get_node_index(column + velocity_vectors[direction][0], VERTICAL_NODES - 1 + velocity_vectors[direction][1]), direction)] = source[access_function(access::get_node_index(column, VERTICAL_NODES - 1), direction)];
        }
    }
}

void stream::two_lattice::helper::two_lattice_wallup(
    all_distributions &source,
    all_distributions &destination, 
    access_function &access_function
    )
{
    for (auto column = 1; column < HORIZONTAL_NODES; ++column)
    {
        for (auto direction : boundaries::neighbor_directions[boundaries::boundary_scenarios::wall_up])
        {
            destination[access_function(access::get_node_index(column + velocity_vectors[direction][0], velocity_vectors[direction][1]), direction)] = source[access_function(access::get_node_index(column, 0), direction)];
        }
    }
}

void stream::two_lattice::helper::two_lattice_regular(
    all_distributions &source,
    all_distributions &destination, 
    access_function &access_function
    )
{
    for (auto row = 1; row < VERTICAL_NODES - 2; ++row)
    {
        for (auto column = 1; column < HORIZONTAL_NODES - 2; ++column)
        {
            for (auto direction : general_stream_directions)
            {
                destination[access_function(access::get_node_index(column + velocity_vectors[direction][0], row + velocity_vectors[direction][1]), direction)] = source[access_function(access::get_node_index(column, row), direction)];
            }
        }
    }
}

void stream::two_lattice::helper perform_boundary_update(all_distributions &source, all_distributions &destination, access_function &access_function)
        {

            // Upper inlet
            boundaries::upper_inlet_boundary_stream(source, destination, access_function);
            // Lower inlet
            boundaries::lower_inlet_boundary_stream(source, destination, access_function);
            // Upper outlet
            boundaries::upper_outlet_boundary_stream(source, destination, access_function);
            // Lower outlet
            boundaries::lower_outlet_boundary_stream(source, destination, access_function);
            // Upper wall and lower wall
            for (int x = 0; x < HORIZONTAL_NODES - 1; ++x)
            {
                boundaries::upper_wall_boundary_stream(source, destination, access_function, x);
                boundaries::lower_wall_boundary_stream(source, destination, access_function, x);
            }
            // inlet and outlet
            for (int y = 0; y < VERTICAL_NODES - 1; ++y)
            {
                boundaries::inlet_boundary_stream(source, destination, access_function, y);
                boundaries::inlet_boundary_stream(source, destination, access_function, y);
            }
        }