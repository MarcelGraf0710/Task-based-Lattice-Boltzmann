#include <iostream>
#include "include/access.hpp"
#include "include/simulation.hpp"
#include "include/utils.hpp"
#include "include/new_two_lattice.hpp"
#include "include/two_step_sequential.hpp"


int main()
{
    // std::vector<double> distribution_values_0;
    // std::vector<unsigned int> nodes;
    // std::vector<unsigned int> fluid_nodes;
    // std::vector<bool> phase_information;
    // border_swap_information swap_info;
    // access_function access_function = access::collision;

    // std::cout << "All vectors declared " << std::endl;

    // setup_example_domain(distribution_values_0, nodes, fluid_nodes, phase_information, swap_info, access_function);

    // std::cout << "distribution_values has size "<< distribution_values_0.size() << std::endl;
    // std::cout << "nodes has size "<< nodes.size() << std::endl;
    // std::cout << "fluid_nodes has size "<< fluid_nodes.size() << std::endl;
    // std::cout << "phase_information has size "<< phase_information.size() << std::endl;

    // std::vector<double> distribution_values_1 = distribution_values_0;

    // std::vector<sim_data_tuple> results = 
    // two_lattice_sequential::run(
    //     fluid_nodes, 
    //     swap_info, 
    //     distribution_values_0, 
    //     distribution_values_1,   
    //     access_function,
    //     5
    //     );

    //    return 0;

    std::vector<double> distribution_values;
    std::vector<unsigned int> nodes;
    std::vector<unsigned int> fluid_nodes;
    std::vector<bool> phase_information;
    border_adjacency ba;
    access_function access_function = access::collision;

    std::cout << "All vectors declared " << std::endl;

    setup_example_domain(distribution_values, nodes, fluid_nodes, phase_information, ba, access_function);

        std::cout << "'distribution_values' has size "<< distribution_values.size() << std::endl;
    std::cout << std::endl;

    std::cout << "Illustration of lattice: " << std::endl;
    to_console::print_phase_vector(phase_information);
    std::cout << std::endl;

    std::cout << "Enumeration of all nodes within the lattice: " << std::endl;
    to_console::print_vector(nodes);
    std::cout << std::endl;

    std::cout << "Enumeration of all fluid nodes within the simulation domain: " << std::endl;
    to_console::print_vector(fluid_nodes, HORIZONTAL_NODES - 2);
    std::cout << std::endl;

    std::cout << "Illustration of border adjacencies: " << std::endl;
    to_console::print_border_adjacencies(ba);
    std::cout << std::endl;

    std::vector<sim_data_tuple> result = two_step_sequential::run(fluid_nodes, distribution_values, ba, access_function, 5);

    return 0;
}