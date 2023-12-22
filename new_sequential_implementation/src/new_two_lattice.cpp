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
sim_data_tuple two_lattice_sequential::perform_tl_stream_and_collide
(
    std::vector<unsigned int> &fluid_nodes,
    border_swap_information &bsi,
    std::vector<double> &source, 
    std::vector<double> &destination,    
    access_function access_function
)
{
    std::set<unsigned int> remaining_nodes = {fluid_nodes.begin(), fluid_nodes.end()};
    std::set<unsigned int> remaining_dirs = {streaming_directions.begin(), streaming_directions.end()};
    std::vector<velocity> velocities(source.size() / DIRECTION_COUNT, velocity{0,0});
    std::vector<double> densities(source.size() / DIRECTION_COUNT, -1);
    std::vector<double> current_distributions(DIRECTION_COUNT, 0);

    std::cout << "\t TL stream and collide: initializations and declarations" << std::endl;

    /* Boundary node treatment */
    bounce_back::perform_early_boundary_update(bsi, destination, access_function);
    std::cout << "\t Early boundary update performed." << std::endl;
    
    for(auto current_border_info : bsi)
    {
        remaining_dirs = two_lattice_sequential::determine_remaining_directions(current_border_info);

        for(auto dir : remaining_dirs) 
            two_lattice_sequential::tl_stream(source, destination, access_function, current_border_info[0], dir);

        current_distributions = access::get_all_distribution_values(destination, current_border_info[0], access_function);
        velocities[current_border_info[0]] = macroscopic::flow_velocity(current_distributions);
        densities[current_border_info[0]] = macroscopic::density(current_distributions);

        two_lattice_sequential::tl_collision(destination, current_border_info[0], access_function, velocities, densities);

        remaining_nodes.erase(current_border_info[0]);
    }

    std::cout << "\t Performed stream and collision for all border nodes." << std::endl;

    /* Treatment of non-boundary nodes */
    for(auto fluid_node : remaining_nodes)
    {   
        for (auto direction = 0; direction < DIRECTION_COUNT; ++direction)
            two_lattice_sequential::tl_stream(source, destination, access_function, fluid_node, direction);

        velocities[fluid_node] = macroscopic::flow_velocity(current_distributions);
        densities[fluid_node] = macroscopic::density(current_distributions);
        two_lattice_sequential::tl_collision(destination, fluid_node, access_function, velocities, densities);
    }
    std::cout << "\t Done treating all non-boundary nodes." << std::endl;
    sim_data_tuple result{velocities, densities};
    return result;
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
void two_lattice_sequential::OLD_tl_collision
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
    std::vector<double >vals = collision::collide_bgk(vals, velocities[fluid_node], densities[fluid_node]);
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
    std::vector<sim_data_tuple> &data
)
{
    std::cout << "--------------------------------------------------------------------------------" << std::endl;
    std::cout << "Now running sequential two lattice algorithm for " << iterations << " iterations." << std::endl;
    std::vector<double> &source = values_0;
    std::vector<double> &destination = values_1;
    std::vector<double> &temp = values_1;
    std::vector<sim_data_tuple>result(iterations, std::make_tuple(std::vector<velocity>(TOTAL_NODE_COUNT, {0,0}), std::vector<double>(TOTAL_NODE_COUNT, 0)));

    for(auto time = 0; time < iterations; ++time)
    {
        std::cout << "Iteration " << time << ":" << std::endl;
        result[time] = two_lattice_sequential::perform_tl_stream_and_collide
        (
            fluid_nodes, 
            boundary_nodes, 
            source, 
            destination, 
            access_function
        );
            
        std::cout << "\tBoth stream and collide performed for all nodes, back at run(...), now changing source and destination..." << std::endl;
        temp = source;
        source = destination;
        destination = temp;
    }
    std::cout << "All done. " << std::endl;
    std::cout << std::endl;


    std::cout << "All done. " << std::endl;
    std::cout << std::endl;
    std::cout << "Velocity values: " << std::endl;
    std::cout << std::endl;
    for(auto i = 0; i < iterations; ++i)
    {
        std::cout << "t = " << i << std::endl;
        std::cout << "-------------------------------------------------------------------------------- " << std::endl;
        print_velocity_vector(std::get<0>(result[i]));
        std::cout << std::endl;
    }
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "Density values: " << std::endl;
    std::cout << std::endl;
    
    for(auto i = 0; i < iterations; ++i)
    {
        std::cout << "t = " << i << std::endl;
        std::cout << "-------------------------------------------------------------------------------- " << std::endl;
        print_vector(std::get<1>(result[i]));
        std::cout << std::endl;
    }
    
}


