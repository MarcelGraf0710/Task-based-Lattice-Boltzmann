#include <iostream>
#include "include/access.hpp"
#include "include/simulation.hpp"
#include "include/utils.hpp"
#include "include/new_two_lattice.hpp"


int main()
{
    std::vector<double> distribution_values_0;
    std::vector<unsigned int> nodes;
    std::vector<unsigned int> fluid_nodes;
    std::vector<bool> phase_information;
    border_swap_information swap_info;
    access_function access_function = access::collision;

    std::cout << "All vectors declared " << std::endl;

    setup_example_domain(distribution_values_0, nodes, fluid_nodes, phase_information, swap_info, access_function);

    std::cout << "distribution_values has size "<< distribution_values_0.size() << std::endl;
    std::cout << "nodes has size "<< nodes.size() << std::endl;
    std::cout << "fluid_nodes has size "<< fluid_nodes.size() << std::endl;
    std::cout << "phase_information has size "<< phase_information.size() << std::endl;

    std::vector<double> distribution_values_1 = distribution_values_0;

    std::vector<sim_data_tuple> results = 
    two_lattice_sequential::run(
        fluid_nodes, 
        swap_info, 
        distribution_values_0, 
        distribution_values_1,   
        access_function,
        5
        );

    return 0;
}