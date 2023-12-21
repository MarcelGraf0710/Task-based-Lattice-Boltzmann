#include "include/access.hpp"
#include "include/simulation.hpp"
#include "include/utils.hpp"
#include "include/new_two_lattice.hpp"
#include "include/two_step_sequential.hpp"

int main()
{
    std::vector<double> distribution_values_0;
    std::vector<unsigned int> nodes;
    std::vector<unsigned int> fluid_nodes;
    std::vector<bool> phase_information;
    sim_data_tuple data{{}, {}};
    border_swap_information swap_info;
    access_function access_function = access::collision;

    std::cout << "All vectors declared. " << std::endl;
    std::cout << std::endl;

    setup_example_domain(distribution_values_0, nodes, fluid_nodes, phase_information, swap_info, access_function);

    std::cout << "'distribution_values' has size "<< distribution_values_0.size() << std::endl;
    std::cout << std::endl;

    std::cout << "Illustration of lattice: " << std::endl;
    print_phase_vector(phase_information);
    std::cout << std::endl;

    std::cout << "Enumeration of all nodes within the lattice: " << std::endl;
    print_vector(nodes);
    std::cout << std::endl;

    std::cout << "Enumeration of all fluid nodes within the simulation domain: " << std::endl;
    print_vector(fluid_nodes, HORIZONTAL_NODES - 2);
    std::cout << std::endl;

    std::cout << "Swap info:" << std::endl;
    for(auto current : swap_info)
        print_vector(current);
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
        data
    );

    return 0;
}