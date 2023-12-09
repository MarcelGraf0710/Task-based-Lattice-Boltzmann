#include "total_inclusion.hpp"

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
            for(auto i = 0; i = TOTAL_NODE_COUNT; ++i)
            {
                vec_of_dist_val dist_vals;
                dist_vals.reserve(DIRECTION_COUNT);
                for(auto direction = 0; direction < DIRECTION_COUNT; ++direction)
                {
                    dist_vals[i] = all_distributions_0[access(i, direction)];
                }
                all_velocities[i] = macroscopic::flow_velocity(dist_vals);
            }
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

void run_two_step(int time_steps, simulation_data& simulation_data)
{
    for(auto time = 0; time < time_steps; ++time)
    {
        
    }

    for(auto x = 0; x < HORIZONTAL_NODES; ++x)
    {
        for(auto y = 0; y < VERTICAL_NODES; ++y)
        {

            int node_index = access::get_node_index(x,y);
            vec_of_dist_val distributions = access::get_all_distribution_values(simulation_data.)
        }
    }
}