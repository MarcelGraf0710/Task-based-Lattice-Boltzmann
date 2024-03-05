#include <iostream>
#include "include/access.hpp"
#include "include/simulation.hpp"
#include "include/utils.hpp"
#include "include/shift_sequential_new.hpp"

int main(const int argc, const char** argv)
{
    bool enable_debug = false;
    if(argc > 2)
    {
        std::cout << "Illegal argument count. Choose argument 'debug' if you want to run this program in debug mode, or run it without arguments." << std::endl;
        exit(-1);
    }
    else if(argc == 2)
    {
        std::string selection(argv[1]);
        enable_debug = debug_handler(selection);
    }

    std::cout << std::endl;
    std::cout << "Starting simulation..." << std::endl;
    std::cout << std::endl;
    std::cout << std::setprecision(3) << std::fixed;

    /* Initializations */
    std::vector<double> distribution_values(0, (TOTAL_NODE_COUNT + SHIFT_OFFSET) * DIRECTION_COUNT);
    std::vector<unsigned int> nodes(0, TOTAL_NODE_COUNT);
    std::vector<unsigned int> fluid_nodes(0, TOTAL_NODE_COUNT);
    std::vector<bool> phase_information(false, TOTAL_NODE_COUNT);
    border_swap_information swap_info;
    access_function access_function = lbm_access::collision;
    
    if(enable_debug)
    {
        std::cout << "All vectors declared. " << std::endl;
        std::cout << std::endl;
    }

    /* Setting up example domain */
    shift_sequential::setup_example_domain(distribution_values, nodes, fluid_nodes, phase_information, swap_info, access_function, enable_debug);

    if(enable_debug)
    {
        /* Illustration of the phase information */
        std::cout << "Illustration of lattice: " << std::endl;
        to_console::print_phase_vector(phase_information);
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
        to_console::print_distribution_values(distribution_values, access_function);
        std::cout << std::endl;
    }

    /* Run simulation */
    shift_sequential::run(fluid_nodes, distribution_values, swap_info, access_function, TIME_STEPS);
}