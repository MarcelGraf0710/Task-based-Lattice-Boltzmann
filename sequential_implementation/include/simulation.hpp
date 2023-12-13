#ifndef SIMULATION_HPP
#define SIMULATION_HPP
#include "defines.hpp"

struct simulation_data
{
    std::vector<double> all_distributions_0;
    std::vector<double> all_distributions_1;
    std::vector<velocity> all_velocities; 
    access_function access;

    simulation_data(std::vector<double> initial_distributions, access_function access);
};

/**
 * @brief Returns a vector containing the distribution values set up according to the specified access pattern and inlet velocity.
 * 
 * @param access 
 * @param inlet_velocity 
 * @return std::vector<double> 
 */
std::vector<double> setup_distributions(access_function access, double inlet_velocity = INLET_VELOCITY, double inlet_density = INLET_DENSITY);

/**
 * @brief
 *
 * @param time_steps
 * @param simulation_data
 * @param access
 */
void run_two_lattice(int time_steps, simulation_data &simulation_data, access_function access);

#endif