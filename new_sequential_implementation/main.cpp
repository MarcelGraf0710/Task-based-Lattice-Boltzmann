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
    std::cout << initial_distributions[0] << std::endl;
    
    simulation_data data;
    std::cout << "Are you ready for some segmentation faults?" << std::endl;
    data.access = access;
    data.all_distributions_0 = initial_distributions;
    data.all_distributions_1 = initial_distributions;
    std::vector<velocity> all_velocities;
    all_velocities.reserve(TOTAL_NODE_COUNT); 
    //macroscopic::update_all_velocities(data.all_distributions_0, all_velocities, access);
    std::cout << &initial_distributions<< std::endl;
    std::cout << &initial_distributions[0] << std::endl;
    std::cout << &initial_distributions[1] << std::endl;
    std::cout << "leaving constructor " << std::endl;
    std::cout << initial_distributions[data.access(0,0)] << std::endl;
    std::cout << &data << std::endl;
    std::cout << &data.all_distributions_0 << std::endl;
    std::cout << &data.all_distributions_1 << std::endl;
    std::cout << &data.access << std::endl;
    std::cout << &data.all_distributions_0[0] << std::endl;
    std::cout << &data.all_distributions_1[0] << std::endl;
    std::cout << "Are you ready for some segmentation faults?" << std::endl;
    std::cout << "Go away" << std::endl;
    std::cout << initial_distributions[0] << std::endl;

    std::cout << data.all_distributions_0[0] << std::endl;
    std::cout << data.all_distributions_1[0] << std::endl;
    std::cout << "Successful access of constructor of simulation data" << std::endl;
    data.all_distributions_1 = data.all_distributions_0;
    std::cout << "Successful setup of simulation data" << std::endl;

    run_two_lattice(5, data, access);
    return 0;
}