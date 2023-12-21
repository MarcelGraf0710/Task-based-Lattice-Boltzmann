#include "../include/new_two_lattice.hpp"
#include "../include/access.hpp"
#include "../include/defines.hpp"
#include "../include/boundaries.hpp"
#include "../include/new_collision.hpp"
#include "../include/macroscopic.hpp"
#include "../include/utils.hpp"
#include <set>
#include <iostream>

/**
 * @brief 
 * 
 * @param fluid_nodes 
 * @param bsi 
 * @param source 
 * @param destination 
 * @param access_function 
 * @param sim_data 
 */
void two_lattice_sequential::perform_tl_stream_and_collide
(
    std::vector<unsigned int> &fluid_nodes,
    border_swap_information &bsi,
    std::vector<double> &source, 
    std::vector<double> &destination,    
    access_function access_function,
    std::vector<double> &sim_data
)
{
    std::set<unsigned int> remaining_nodes = {fluid_nodes.begin(), fluid_nodes.end()};
    std::set<unsigned int> remaining_dirs;
    std::vector<velocity> velocities(fluid_nodes.size(), velocity{0,0});
    std::vector<double> densities(fluid_nodes.size(), 0);
    std::cout << "\t TL stream and collide: initializations and declarations" << std::endl;

    /* Boundary node treatment */
    bounce_back::perform_early_boundary_update(bsi, destination, access_function);
    std::cout << "\t Early boundary update performed." << std::endl;
    
    for(auto current_border_info : bsi)
    {
        std::cout << "Currently dealing with node " << current_border_info[0] << std::endl;
        two_lattice_sequential::determine_remaining_directions(current_border_info, remaining_dirs);
        for(auto dir : remaining_dirs) 
            two_lattice_sequential::tl_stream(source, destination, access_function, current_border_info[0], dir);
        two_lattice_sequential::tl_collision(destination, current_border_info[0], access_function, velocities, densities);
        remaining_nodes.erase(current_border_info[0]);
    }

    std::cout << "\t Performed stream and collision for all border nodes." << std::endl;

    /* Treatment of non-boundary nodes */
    for(auto fluid_node : remaining_nodes)
    {   
        for (auto direction = 0; direction < DIRECTION_COUNT; ++direction)
            two_lattice_sequential::tl_stream(source, destination, access_function, fluid_node, direction);
        two_lattice_sequential::tl_collision(destination, fluid_node, access_function, velocities, densities);
    }
    std::cout << "Done treating all non-boundary nodes." << std::endl;

    std::vector<double> new_sim_data(10,0); //get_simulation_data(velocities, densities);
    sim_data.insert(sim_data.end(), new_sim_data.begin(), new_sim_data.end());
    //std::cout << "Refernces set." << std::endl;

    //print_vector(velocities_sim_data);
    //print_vector(densities, HORIZONTAL_NODES - 2);
    //print_velocity_vector(velocities);

    //std::vector<velocity> test = std::move(velocities);
    //std::cout << "First: (" << (*velocities.begin())[0] << ", " << (*velocities.begin())[1] << ")" << std::endl;
    //std::cout << "Last: (" << (*velocities.end()--)[0] << ", " << (*velocities.end()--)[1] << ")" << std::endl;
    //print_velocity_vector(test);
    //std::cout << *velocities.begin() << std::endl;
    //velocities_sim_data.insert(velocities_sim_data.begin(), velocities.begin(), velocities.end());
    std::cout << "Velocities set." << std::endl;
    std::cout << "Densities set." << std::endl;
}

/**
 * @brief 
 * 
 * @param destination 
 * @param fluid_node 
 * @param access_function 
 * @param velocities 
 * @param densities 
 */
void two_lattice_sequential::tl_collision
(
    std::vector<double> &destination, 
    unsigned int fluid_node, 
    access_function &access_function, 
    std::vector<velocity> &velocities, 
    std::vector<double> &densities
)
{
    std::vector<double> vals = access::get_all_distribution_values(destination, fluid_node, access_function);
    velocity current_velocity = macroscopic::flow_velocity(vals);
    velocities[fluid_node] = current_velocity;
    double density = macroscopic::density(vals);
    densities[fluid_node] = density;
    vals = collision::collide_bgk(vals, current_velocity, density);
    access::set_all_distribution_values(vals, destination, fluid_node, access_function);
}


void two_lattice_sequential::run
(  
    std::vector<unsigned int> &fluid_nodes,       
    border_swap_information &boundary_nodes,
    std::vector<double> &values_0, 
    std::vector<double> &values_1,   
    access_function access_function,
    unsigned int iterations,
    std::vector<double> &data
)
{
    std::cout << "--------------------------------------------------------------------------------" << std::endl;
    std::cout << "Now running sequential two lattice algorithm for " << iterations << " iterations." << std::endl;
    std::vector<double> &source = values_0;
    std::vector<double> &destination = values_1;
    std::vector<double> &temp = values_1;
    std::vector<std::vector<double>> result;

    for(auto time = 0; time < iterations; ++time)
    {
        std::cout << "Iteration " << time << ":" << std::endl;
        two_lattice_sequential::perform_tl_stream_and_collide
        (
            fluid_nodes, 
            boundary_nodes, 
            source, 
            destination, 
            access_function, 
            data
        );
            
        std::cout << "\tBoth stream and collide performed for all nodes, back at run(...), now changing source and destination..." << std::endl;
        temp = source;
        source = destination;
        destination = temp;
    }
    std::cout << "All done. " << std::endl;
    
}


