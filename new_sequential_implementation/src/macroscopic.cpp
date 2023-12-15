#include "../include/macroscopic.hpp"
#include <iostream>

double macroscopic::density(vec_of_dist_val &distribution_functions)
{
    return std::accumulate(distribution_functions.begin(), distribution_functions.end(), 0);
}

velocity macroscopic::flow_velocity(vec_of_dist_val &distribution_functions)
{
    double v_x, v_y;
    velocity velocity_vector;

    for(int i = 0; i < DIRECTION_COUNT; ++i)
    {
        velocity_vector = velocity_vectors[i];
        v_x += distribution_functions[i] * velocity_vector[0];
        v_y += distribution_functions[i] * velocity_vector[1];
    }
    std::cout << "Leaving flow velocity " << std::endl;
    velocity v{v_x, v_y};
    return v;
}

void macroscopic::update_all_velocities(std::vector<double> &all_distributions, std::vector<velocity> &destination, access_function access_f)
{
    for(auto x = 0; x < HORIZONTAL_NODES; ++x)
    {
        for(auto y = 0; y < VERTICAL_NODES; ++y)
        {
            int node_index = access::get_node_index(x,y);
            std::vector<double> dist_vals = access::get_all_distribution_values(all_distributions, node_index, access_f);
            std::cout << "I'm back " << std::endl;
            std::cout << dist_vals[0] << std::endl;
            std::cout << "Flow velocity is " << flow_velocity(dist_vals)[0] << flow_velocity(dist_vals)[1] << std::endl;
            std::cout << "(" << destination[node_index][0] << ", " << destination[node_index][1] << ")" << std::endl;
            destination[node_index] = flow_velocity(dist_vals);
            std::cout << "Setup destination value of (" << x << ", " << y << ")" << std::endl;
        }
    }
    std::cout << "Updated all velocities " << std::endl;
}

