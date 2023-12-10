#include "defines.hpp"
#include "macroscopic.hpp"
#include "access.hpp"
#include "update.hpp"

struct simulation_data
{
    std::vector<double> all_distributions_0;
    std::vector<double> all_distributions_1;
    std::vector<velocity> all_velocities; 
    access_function access;

    simulation_data(std::vector<double> initial_distributions, access_function access) : 
        all_distributions_0(initial_distributions),
        access(access)
        {
            macroscopic::update_all_velocities(all_distributions_0, all_velocities, access);
        };
};

/**
 * @brief Returns a vector containing the distribution values set up according to the specified access pattern and inlet velocity.
 * 
 * @param access 
 * @param inlet_velocity 
 * @return std::vector<double> 
 */
std::vector<double> setup_distributions(access_function access, double inlet_velocity = INLET_VELOCITY, double inlet_density = INLET_DENSITY)
{
    // Setup all non-inlet nodes to have initial velocity 0
    std::vector<double> dist_vals;
    dist_vals.reserve(TOTAL_NODE_COUNT * DIRECTION_COUNT);
    for(auto x = 1; x < HORIZONTAL_NODES; ++x)
    {
        for(auto y = 0; y < VERTICAL_NODES; ++y)
        {
            for(auto direction = 0; direction < DIRECTION_COUNT; ++direction)
            {
                dist_vals[access(access::get_node_index(x,y), direction)] = maxwell_boltzmann_distribution({0,0}, inlet_density, direction);
            }
        }
    }
    // Setup all inlet nodes to have initial velocity 0
    for(auto y = 0; y < VERTICAL_NODES; ++y)
        {
            for(auto direction = 0; direction < DIRECTION_COUNT; ++direction)
            {
                dist_vals[access(access::get_node_index(0,y), direction)] = maxwell_boltzmann_distribution({inlet_velocity,0}, inlet_density, direction);
            }
        }

    return dist_vals; 
}

void run_two_lattice(int time_steps, simulation_data& simulation_data, access_function access)
{
    for(auto time = 0; time < time_steps; ++time)
    {
        all_distributions& source = (time % 2) ? simulation_data.all_distributions_0 : simulation_data.all_distributions_1;
        all_distributions& destination = (time % 2) ? simulation_data.all_distributions_1 : simulation_data.all_distributions_0;
        collision::perform_collision_step(source, simulation_data.all_velocities, access);
        stream::two_lattice::perform_two_lattice_stream(access, source, destination);
        macroscopic::update_all_velocities(destination, simulation_data.all_velocities, access);
    }
}