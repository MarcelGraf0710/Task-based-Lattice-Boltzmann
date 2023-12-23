#include "../include/update.hpp"
#include <iostream>

std::array<std::tuple<std::tuple<unsigned int, unsigned int>, std::list<int>>, 4> stream::corner_specs{
    std::make_tuple(std::make_tuple(0,0), std::list{5,7,8}),
    std::make_tuple(std::make_tuple(0,VERTICAL_NODES - 1), std::list{1,2,5}),
    std::make_tuple(std::make_tuple(HORIZONTAL_NODES - 1,0), std::list{3,6,7}),
    std::make_tuple(std::make_tuple(HORIZONTAL_NODES - 1,VERTICAL_NODES - 1), std::list{0,1,3})
};

std::array<unsigned int, DIRECTION_COUNT - 1> stream::general_stream_directions{0,1,2,3,5,6,7,8};

vec_of_dist_val collision::collide_bgk(vec_of_dist_val &f, velocity u, double density)
{
    std::cout << "Entering collision" << std::endl;
    vec_of_dist_val result = maxwell_boltzmann_distribution(u, density);
    std::cout << "Got mb distri" << std::endl;
    for(auto i = 0; i < DIRECTION_COUNT; ++i)
    {
        std::cout << "Loop iteration " << i << std::endl;
        result[i] = -(1/RELAXATION_TIME) * (f[i] - result[i]);
    }
    return result;
}

void collision::perform_collision_step(
    std::vector<double> &all_distribution_values, 
    std::vector<velocity> &all_velocities, 
    access_function access, 
    double density)
{
    for(auto x = 0; x < HORIZONTAL_NODES; ++x)
    {
        for(auto y = 0; y < VERTICAL_NODES; ++y)
        {
            int node_index = access::get_node_index(x,y);
            std::cout << all_distribution_values[0] << std::endl;
            std::vector<double> current_dist_values = access::get_distribution_values_of(all_distribution_values, node_index, access);
            std::cout << "Set current dist values, entering collide_bgk" << std::endl;
            current_dist_values = collide_bgk(current_dist_values, all_velocities[node_index], density);
            std::cout << "Setting all distribution values" << std::endl;
            access::set_all_distribution_values(current_dist_values, all_distribution_values, node_index, access);
        }
    }
}

void stream::two_lattice::helper::two_lattice_corners(
    all_distributions &source,
    all_distributions &destination, 
    access_function &access_function
    )
{
    for (auto corner : stream::corner_specs)
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
        for (auto direction : std::list{0,1,3,6,7})
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
        for (auto direction : std::list{1,2,5,7,8})
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
        for (auto direction : std::list{3,5,6,7,8})
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
        for (auto direction : std::list{0,1,2,3,5})
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
            for (auto direction : stream::general_stream_directions)
            {
                destination[access_function(access::get_node_index(column + velocity_vectors[direction][0], row + velocity_vectors[direction][1]), direction)] = source[access_function(access::get_node_index(column, row), direction)];
            }
        }
    }
}

void stream::two_lattice::helper::perform_boundary_update(all_distributions &source, all_distributions &destination, access_function &access_function)
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

void stream::two_lattice::perform_two_lattice_stream(
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