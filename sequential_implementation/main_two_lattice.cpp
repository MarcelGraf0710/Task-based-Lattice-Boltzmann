#include "include/access.hpp"
#include "include/simulation.hpp"
#include "include/utils.hpp"
#include "include/two_lattice_sequential.hpp"

int main()
{
    std::cout << "Starting simulation..." << std::endl;

    /* Initializations */
    std::vector<double> distribution_values_0(0, TOTAL_NODE_COUNT * DIRECTION_COUNT);
    std::vector<unsigned int> nodes(0, TOTAL_NODE_COUNT);
    std::vector<unsigned int> fluid_nodes(0, TOTAL_NODE_COUNT);
    std::vector<bool> phase_information(false, TOTAL_NODE_COUNT);
    std::vector<sim_data_tuple>result(5, std::make_tuple(std::vector<velocity>(TOTAL_NODE_COUNT, {0,0}), std::vector<double>(TOTAL_NODE_COUNT, 0)));
    border_swap_information swap_info;
    access_function access_function = access::bundle;
    std::cout << "All vectors declared. " << std::endl;
    std::cout << std::endl;

    /* Setting up example domain */
    setup_example_domain(distribution_values_0, nodes, fluid_nodes, phase_information, swap_info, access_function);

    /* Illustration of the phase information */
    std::cout << "Illustration of lattice: " << std::endl;
    to_console::print_phase_vector(phase_information);
    std::cout << std::endl;

    /* Color information */
    std::cout << "This program utilizes ANSI color codes to output colored text. If your command line does not support those codes, your output may be corrupted." << std::endl;
    std::cout << "In all following prints showing the entire simulation domain, "; 
    std::cout << "the origin will be marked in \033[31mred\033[0m and the outmost coordinate will be marked in \033[34mblue\033[0m." << std::endl;
    std::cout << "Milestones will be marked in \033[33myellow\033[0m." << std::endl;
    std::cout << std::endl;

    /* Overview */
    std::cout << "Enumeration of all nodes within the lattice: " << std::endl;
    to_console::print_vector(nodes);
    std::cout << std::endl;

    std::cout << "Enumeration of all fluid nodes within the simulation domain: " << std::endl;
    to_console::print_vector(fluid_nodes, HORIZONTAL_NODES - 2);
    std::cout << std::endl;

    std::cout << "Swap info:" << std::endl;
    for(const auto& current : swap_info)
        to_console::print_vector(current, current.size());
    std::cout << std::endl;

    std::cout << "Initial distributions:" << std::endl;
    to_console::print_distribution_values(distribution_values_0, access_function);
    std::cout << std::endl;

    std::vector<double> distribution_values_1 = distribution_values_0;

    /* Run simulation */
    two_lattice_sequential::run
    (
        fluid_nodes, 
        swap_info, 
        distribution_values_0, 
        distribution_values_1,   
        access_function,
        10,
        result
    );

    return 0;
}