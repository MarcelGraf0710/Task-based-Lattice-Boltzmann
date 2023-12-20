#include "../include/new_two_lattice.hpp"
#include "../include/access.hpp"
#include "../include/defines.hpp"
#include "../include/boundaries.hpp"
#include "../include/new_collision.hpp"
#include "../include/macroscopic.hpp"
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
    sim_data_vector &sim_data
)
{
    std::set<unsigned int> remaining_nodes = {fluid_nodes.begin(), fluid_nodes.end()};
    std::set<unsigned int> remaining_dirs;
    std::vector<velocity> velocities(fluid_nodes.size(), velocity{0,0});
    std::vector<double> densities(fluid_nodes.size(), 0);

    /* Boundary node treatment */
    bounce_back::perform_early_boundary_update(bsi, destination, access_function);
    std::cout << "\t Early boundary update performed." << std::endl;

    for(auto current_border_info : bsi)
    {
        remaining_dirs = two_lattice_sequential::determine_remaining_directions(current_border_info);
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

    std::vector<velocity> &velocities_sim_data = sim_data[0];
    std::vector<double> densities_sim_data = sim_data[1];

    velocities_sim_data.insert(velocities_sim_data.end(), velocities.begin(), velocities.end());
    densities_sim_data.insert(densities_sim_data.end(), densities.begin(), densities.end());
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
    sim_data_vector &data
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


