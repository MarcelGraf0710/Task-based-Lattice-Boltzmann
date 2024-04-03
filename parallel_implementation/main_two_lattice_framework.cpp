#include "include/access.hpp"
#include "include/simulation.hpp"
#include "include/utils.hpp"
#include "include/two_lattice_sequential.hpp"
#include "include/parallel_two_lattice_framework.hpp"
#include "include/parallel_framework.hpp"
#include <hpx/hpx_init.hpp>

int hpx_main()
{
    std::cout << std::endl;
    std::cout << "Starting simulation..." << std::endl;
    std::cout << std::endl;
    std::cout << std::setprecision(3) << std::fixed;

    /* Color message */
    to_console::print_ansi_color_message();


    /* Initializations */
    std::vector<double> distribution_values_0(0, TOTAL_NODE_COUNT * DIRECTION_COUNT);
    std::vector<unsigned int> nodes(0, TOTAL_NODE_COUNT);
    std::vector<unsigned int> fluid_nodes(0, TOTAL_NODE_COUNT);
    std::vector<bool> phase_information(false, TOTAL_NODE_COUNT);
    border_swap_information swap_info;
    access_function access_function = lbm_access::collision;

    std::cout << "All vectors declared. " << std::endl;
    std::cout << std::endl;

    /* Setting up example domain */
    parallel_framework::setup_parallel_domain(distribution_values_0, nodes, fluid_nodes, phase_information, access_function);

    std::vector<start_end_it_tuple> subdomain_fluid_bounds;
    for(auto subdomain = 0; subdomain < SUBDOMAIN_COUNT; ++subdomain)
    {
        // std::cout << "Currently determining bounds for subdomain " << subdomain << std::endl;
        subdomain_fluid_bounds.push_back(parallel_framework::get_subdomain_fluid_node_pointers(subdomain, fluid_nodes));
    }

    swap_info = parallel_framework::retrieve_fast_border_swap_info(subdomain_fluid_bounds, fluid_nodes, phase_information);

    // /* Illustration of the phase information */
    // std::cout << "Illustration of lattice: " << std::endl;
    // to_console::print_phase_vector(phase_information);
    // std::cout << std::endl;

    // /* Overview */
    // std::cout << "Enumeration of all nodes within the lattice: " << std::endl;
    // to_console::print_vector_buffered(nodes);
    // std::cout << std::endl;

    // std::cout << "Enumeration of all fluid nodes within the simulation domain: " << std::endl;
    // to_console::print_vector(fluid_nodes, HORIZONTAL_NODES - 2);
    // std::cout << std::endl;

    // std::cout << "Swap info:" << std::endl;
    // for(const auto& current : swap_info)
    //     to_console::print_vector(current, current.size());
    // std::cout << std::endl;

    // std::cout << "Initial distributions:" << std::endl;
    // to_console::print_distribution_values_buffered(distribution_values_0, access_function);
    // std::cout << std::endl;

    std::vector<double> distribution_values_1 = distribution_values_0;

    /* Run simulation */
    parallel_two_lattice_framework::run
    (
        subdomain_fluid_bounds, 
        swap_info, 
        distribution_values_0, 
        distribution_values_1,   
        access_function,
        TIME_STEPS
    );

    return hpx::finalize();
}

int main()
{
    // Initialize HPX, run hpx_main as the first HPX thread, and
    // wait for hpx::finalize being called.
    return hpx::init();
}