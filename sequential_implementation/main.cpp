#include <iostream>
#include "access.hpp"
#include "boundaries.hpp"
#include "defines.hpp"
#include "macroscopic.hpp"
#include "simulation.cpp"
#include "update.hpp"

int main()
{
    std::cout << "Hello World!\n" << std::endl;
    access_function access = access::bundle;
    std::vector<double> initial_distributions = setup_distributions(access);
    simulation_data simulation_data(initial_distributions, access);
    simulation_data.all_distributions_1 = simulation_data.all_distributions_0;
    run_two_lattice(5, simulation_data, access);
    return 0;
}