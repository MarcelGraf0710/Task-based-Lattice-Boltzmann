#include "include/access.hpp"
#include "include/simulation.hpp"
#include "include/utils.hpp"
#include "include/new_two_lattice.hpp"
#include "include/two_step_sequential.hpp"

int main()
{
    std::cout << "Starting simulation..." << std::endl;

    std::vector<double> distribution_values_0;
    distribution_values_0.reserve(TOTAL_NODE_COUNT * DIRECTION_COUNT);
    std::vector<unsigned int> nodes;
    nodes.reserve(TOTAL_NODE_COUNT);
    std::vector<unsigned int> fluid_nodes;
    fluid_nodes.reserve(TOTAL_NODE_COUNT);
    std::vector<bool> phase_information;
    phase_information.reserve(TOTAL_NODE_COUNT);
    std::vector<sim_data_tuple>result(5, std::make_tuple(std::vector<velocity>(TOTAL_NODE_COUNT, {0,0}), std::vector<double>(TOTAL_NODE_COUNT, 0)));
    border_swap_information swap_info;
    access_function access_function = access::collision;

    std::cout << std::endl;
    std::cout << "All vectors declared. " << std::endl;

    setup_example_domain(distribution_values_0, nodes, fluid_nodes, phase_information, swap_info, access_function);

    std::cout << "'distribution_values' has size "<< distribution_values_0.size() << std::endl;
    std::cout << std::endl;

    std::cout << "Illustration of lattice: " << std::endl;
    to_console::print_phase_vector(phase_information);
    std::cout << std::endl;

    std::cout << "This program utilizes ANSI color codes to output colored text. If your command line does not support those codes, your output may be corrupted." << std::endl;
    std::cout << "In all following prints showing the entire simulation domain, "; 
    std::cout << "the origin will be marked in \033[31mred\033[0m and the outmost coordinate will be marked in \033[34mblue\033[0m." << std::endl;
    std::cout << "Milestones will be marked in \033[33myellow\033[0m." << std::endl;
    std::cout << std::endl;

    std::cout << "Enumeration of all nodes within the lattice: " << std::endl;
    to_console::print_vector(nodes);
    std::cout << std::endl;

    std::cout << "Enumeration of all fluid nodes within the simulation domain: " << std::endl;
    to_console::print_vector(fluid_nodes, HORIZONTAL_NODES - 2);
    std::cout << std::endl;

    std::cout << "Swap info:" << std::endl;
    for(auto current : swap_info)
        to_console::print_vector(current, current.size());
    std::cout << std::endl;

    std::cout << "Initial distributions:" << std::endl;
    to_console::print_distribution_values(distribution_values_0, access_function);
    std::cout << std::endl;


    std::vector<double> distribution_values_1 = distribution_values_0;

    two_lattice_sequential::run
    (
        fluid_nodes, 
        swap_info, 
        distribution_values_0, 
        distribution_values_1,   
        access_function,
        5,
        result
    );

    return 0;
}