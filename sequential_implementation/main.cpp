#include <iostream>
#include "include/access.hpp"
//#include "include/simulation_in_cool.hpp"
#include "include/simulation.hpp"

int main()
{
    //print_happy_life();
    std::cout << "Hello World!\n" << std::endl;
    access_function access = access::bundle;
    std::cout << "Successful setup of access function" << std::endl;
    std::vector<double> initial_distributions = setup_distributions(access);
    std::cout << "Successful setup of distribution function" << std::endl;
    simulation_data simulation_data(initial_distributions, access);
    std::cout << "Successful access of constructor of simulation data" << std::endl;
    simulation_data.all_distributions_1 = simulation_data.all_distributions_0;
    std::cout << "Successful setup of simulation data" << std::endl;
    run_two_lattice(5, simulation_data, access);
    return 0;
}