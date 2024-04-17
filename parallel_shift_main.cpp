#include "include/parallel_shift_framework.hpp"

#include <hpx/hpx_init.hpp>

int hpx_main()
{
    bool enable_debug = false;

    std::cout << std::endl;
    std::cout << "Starting simulation..." << std::endl;
    std::cout << std::endl;
    std::cout << std::setprecision(3) << std::fixed;

    /* Color message */
    to_console::print_ansi_color_message();


    /* Initializations */
    std::vector<double> distribution_values;
    std::vector<unsigned int> nodes(0, TOTAL_NODE_COUNT);
    std::vector<unsigned int> fluid_nodes(0, TOTAL_NODE_COUNT);
    std::vector<bool> phase_information(false, TOTAL_NODE_COUNT);
    std::vector<border_swap_information> swap_info;
    access_function access_function = parallel_shift_framework::access_functions::bundle;

    std::cout << "All vectors declared. " << std::endl;
    std::cout << std::endl;

    /* Setting up example domain */
    parallel_shift_framework::setup_parallel_domain(distribution_values, nodes, fluid_nodes, phase_information, access_function);

    std::vector<start_end_it_tuple> subdomain_fluid_bounds;
    for(auto subdomain = 0; subdomain < SUBDOMAIN_COUNT; ++subdomain)
    {
        subdomain_fluid_bounds.push_back(parallel_framework::get_subdomain_fluid_node_pointers(subdomain, fluid_nodes));
    }

    swap_info = parallel_framework::subdomain_wise_border_swap_info(subdomain_fluid_bounds, fluid_nodes, phase_information);

    if(enable_debug)
    {
        /* Illustration of the phase information */
        std::cout << "Illustration of lattice: " << std::endl;
        to_console::print_phase_vector(phase_information);
        std::cout << std::endl;

        /* Overview */
        std::cout << "Enumeration of all nodes within the lattice: " << std::endl;
        to_console::buffered::print_vector(nodes);
        std::cout << std::endl;

        std::cout << "Enumeration of all fluid nodes within the simulation domain: " << std::endl;
        to_console::print_vector(fluid_nodes, HORIZONTAL_NODES - 2);
        std::cout << std::endl;

        std::cout << "Subdomain-wise border swap information: " << std::endl;
        for(auto subdomain = 0; subdomain < SUBDOMAIN_COUNT; ++subdomain)
        {
            std::cout << "Subdomain " <<  subdomain << ":" << std::endl;
            for(const auto& current : swap_info[subdomain])
            to_console::print_vector(current, current.size());
        }
        std::cout << std::endl;
    }

    /* Run simulation */
    parallel_shift_framework::run
    (
        subdomain_fluid_bounds, 
        swap_info, 
        distribution_values, 
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