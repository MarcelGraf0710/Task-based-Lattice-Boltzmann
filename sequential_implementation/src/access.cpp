#include "../include/access.hpp"
#include <iostream>

std::vector<double> access::get_all_distribution_values(std::vector<double> &source, int node_index, access_function access)
{
    std::cout << "Getting all distribution values " << std::endl;
    std::vector<double> dist_vals;
    dist_vals.reserve(9);
    for(auto direction = 0; direction < DIRECTION_COUNT; ++direction)
    {
        std::cout << "Attempting dist_vals[" << direction << "] = ";
        std::cout << source[access(node_index, direction)];
        dist_vals[direction] = source[access(node_index, direction)];
        std::cout << dist_vals[direction] << std::endl;
    }
    return dist_vals;
}

void access::set_all_distribution_values(std::vector<double> &dist_vals, std::vector<double> &destination, int node_index, access_function access)
{
    for(auto direction = 0; direction < DIRECTION_COUNT; ++direction)
    {
        std::cout << "Attempting destination[";
        std::cout << access(node_index, direction);
        std::cout << "] = dist_vals[" << direction << "]" << std::endl;
        std::cout << destination[access(node_index, direction)] << std::endl;
        std::cout << dist_vals[0] << std::endl;
        destination[access(node_index, direction)] = dist_vals[direction];
    }
}

std::vector<unsigned int> semi_direct_access::get_fluid_segments(std::vector<bool> node_phases)
{
    unsigned int index = 0;
    unsigned int consecution = 0;
    std::list<unsigned int> fluid_segments;

    while (index < size(node_phases))
    {
        if (!node_phases[index + consecution]) consecution++; // Hit fluid node
        else if (consecution > 0) // Hit solid node that marks the end of a line of consecutive fluid nodes
        {
            fluid_segments.push_back(index);
            fluid_segments.push_back(consecution);
            index += consecution;
            consecution = 0;
        }
        else index++; // Hit consecutive solid nodes
    }

    std::vector<unsigned int> result
        {
            std::make_move_iterator(begin(fluid_segments)), 
            std::make_move_iterator(end(fluid_segments))
        };

    return result;
}
