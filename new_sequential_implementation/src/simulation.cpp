#include "../include/simulation.hpp"
#include "../include/access.hpp"
#include "../include/macroscopic.hpp"
#include "../include/update.hpp"
#include <iostream>

// simulation_data::simulation_data(const std::vector<double> &initial_distributions, access_function access) : 
//     all_distributions_0(initial_distributions),
//     access(access)
//     {
//         all_distributions_0 = initial_distributions;
//         all_distributions_1 = initial_distributions;
//         std::cout << "Accessing constructor " << std::endl;
//         std::cout << "!!!!!" << initial_distributions[0] << std::endl;
//         all_velocities.reserve(TOTAL_NODE_COUNT); 
//         macroscopic::update_all_velocities(all_distributions_0, all_velocities, access);
//         std::cout << "leaving constructor " << std::endl;
//     };

std::vector<double> setup_distributions(access_function access, double inlet_velocity, double inlet_density)
{
    // Setup all non-inlet nodes to have initial velocity 0
    std::vector<double> dist_vals;
    dist_vals.reserve(TOTAL_NODE_COUNT * DIRECTION_COUNT);
    velocity zero_velocity{0,0};
    velocity inlet_velocity_vector = {inlet_velocity,0};
    for(auto x = 1; x < HORIZONTAL_NODES; ++x)
    {
        for(auto y = 0; y < VERTICAL_NODES; ++y)
        {
            for(auto direction = 0; direction < DIRECTION_COUNT; ++direction)
            {
                double arg = maxwell_boltzmann_distribution(zero_velocity, inlet_density, direction);
                dist_vals[access(access::get_node_index(x,y), direction)] = arg;
            }
        }
    }
    // Setup all inlet nodes to have initial velocity 0
    for(auto y = 0; y < VERTICAL_NODES; ++y)
        {
            for(auto direction = 0; direction < DIRECTION_COUNT; ++direction)
            {
                dist_vals[access(access::get_node_index(0,y), direction)] = maxwell_boltzmann_distribution(inlet_velocity_vector, inlet_density, direction);
            }
        }

    return dist_vals; 
}

void run_two_lattice(int time_steps, simulation_data& simulation_data, access_function access)
{
    std::cout << "Accessing run_two_lattice " << std::endl;
    for(auto time = 0; time < time_steps; ++time)
    {
        std::cout << simulation_data.all_distributions_0[0] << std::endl;
        std::cout << simulation_data.all_distributions_1[0] << std::endl;
        std::cout << "Entering time step " << time << std::endl;
        all_distributions source = (time % 2) ? simulation_data.all_distributions_0 : simulation_data.all_distributions_1;
        all_distributions destination = (time % 2) ? simulation_data.all_distributions_1 : simulation_data.all_distributions_0;
        
        std::cout << simulation_data.all_velocities[0][1] << std::endl;
        std::cout << "Performing collision step " << std::endl;

        collision::perform_collision_step(source, simulation_data.all_velocities, access);
        std::cout << "Performing streaming step " << std::endl;
        stream::two_lattice::perform_two_lattice_stream(access, source, destination);
        std::cout << "Updating velocities " << std::endl;
        macroscopic::calculate_all_velocities(destination, simulation_data.all_velocities, access);
    }
}