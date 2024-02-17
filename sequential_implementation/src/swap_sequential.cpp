#include "../include/swap_sequential.hpp"
#include "../include/access.hpp"
#include "../include/boundaries.hpp"
#include "../include/macroscopic.hpp"
#include "../include/collision.hpp"
#include <iostream>

/**
 * @brief Maps the border node types to the swap partner directions.
 *        0 - {}        (upper right corner node)
 *        1 - {5}       (upper edge node or left upper corner node)
 *        2 - {6,7}     (right edge node or right lower corner node)
 *        3 - {5,7,8}   (left edge node or left lower corner node)
 *        4 - {5,6,7,8} (non-border node or lower edge node)
 */
std::map<char, std::set<unsigned int>> swap_sequential::swap_directions
{
    {0, {}},
    {1, {5}},
    {2, {6, 7}},
    {3, {5, 7, 8}},
    {4, {5, 6, 7, 8}}
};

/**
 * @brief This map stores the inlet and outlet directions for each node.
 *        0 - {}        (empty set, neither inlet nor outlet)
 *        1 - {2,5,8}   (inlet node)
 *        2 - {0,3,6}   (outlet node)
 */
std::map<char, std::set<unsigned int>> swap_sequential::inout_directions
{
    {0, {}},
    {1, {2, 5, 8}},
    {2, {0, 3, 6}}
};

/**
 * @brief Sets up the swap information data structure
 * 
 * @param bsi see documentation of border_swap_information
 * @param fluid_nodes a vector containing the indices of all fluid nodes
 * @return see documentation of swap_information 
 */
swap_information swap_sequential::setup_swap_information
(
    const border_swap_information &bsi, 
    const std::vector<unsigned int> &fluid_nodes
)
{
    unsigned int x = 0;
    unsigned int y = 0;
    unsigned int node = 0;
    std::tuple<unsigned int, unsigned int> coords = {0,0};
    unsigned int bsi_index = 0;
    swap_information result;
    result.reserve(fluid_nodes.size());

    for(auto i = 0; i < fluid_nodes.size(); ++i)
    {
        node = fluid_nodes[i];
        if(node == bsi[bsi_index][0]) // node is at border
        {
            bsi_index++;
            coords = access::get_node_coordinates(node);
            x = std::get<0>(coords);
            y = std::get<1>(coords);

            if(x == 1) // inlet node
            {
                if(y == 1) result.push_back({node,3,1}); // Left lower node
                else if (y == (VERTICAL_NODES - 2)) result.push_back({node,1,1}); // Left upper node
                else result.push_back({node,3,1}); // Left edge node
            }
            else if(x == (HORIZONTAL_NODES - 2)) // outlet node
            {
                if(y == 1) result.push_back({node,2,2}); // Right lower node
                else if (y == (VERTICAL_NODES - 2)) result.push_back({node,0,2}); // Right edge node
                else result.push_back({node,2,2}); // Right upper node
            }
            else // neither inlet nor outlet
            {
                if(y == 1) result.push_back({node,4,0}); // Upper edge node
                else result.push_back({node,1,0}); // lower edge node
            }
        }
        else
        {
            result.push_back({node,4,0});
        }
    }
    return result;
}

/**
 * @brief Performs the combined streaming and collision step for all fluid nodes within the simulation domain.
 *        The border conditions are enforced through ghost nodes.
 *        This variant of the combined streaming and collision step will print several debug comments to the console.
 */
sim_data_tuple swap_sequential::perform_swap_stream_and_collide_debug
(
    const swap_information &swap_information,
    std::vector<double> &distribution_values,    
    const access_function access_function
)
{
    std::cout << "Distribution values before stream and collide: " << std::endl;
    to_console::print_distribution_values(distribution_values, access_function);
    std::cout << std::endl;

    std::vector<velocity> velocities(TOTAL_NODE_COUNT, velocity{0,0});
    velocity current_velocity = {0,0};

    std::vector<double> densities(TOTAL_NODE_COUNT, -1);
    double current_density = -1;

    unsigned int node_index = 0;
    std::set<unsigned int> swap_dirs = {0,0,0};
    std::set<unsigned int> inout = {0,0,0};
    std::vector<double> current_distributions(DIRECTION_COUNT, 0);

    // Swapping step
    for(const auto& node : swap_information)
    {
        node_index = std::get<0>(node);
        swap_dirs = swap_directions[std::get<1>(node)];
        inout = inout_directions[std::get<2>(node)];

        swap_sequential::perform_swap_step(distribution_values, node_index, access_function, swap_dirs);
    }
    std::cout << "Distribution values after swap for every node: " << std::endl;
    to_console::print_distribution_values(distribution_values, access_function);

    // Restore precious order in here
    for(const auto& node : swap_information)
    {
        node_index = std::get<0>(node);
        swap_dirs = swap_directions[std::get<1>(node)];
        inout = inout_directions[std::get<2>(node)];

        swap_sequential::restore_order(distribution_values, node_index, access_function);
    }
    std::cout << "Distribution values after ORDER has been restored for every node: " << std::endl;
    to_console::print_distribution_values(distribution_values, access_function);

    // Perform inflow and outflow, if necessary
    for(const auto& node : swap_information)
    {
        node_index = std::get<0>(node);
        swap_dirs = swap_directions[std::get<1>(node)];
        inout = inout_directions[std::get<2>(node)];

        boundary_conditions::single_node_inout(distribution_values, node_index, inout, access_function);
    }
    std::cout << "Distribution values after inflow and outflow update: " << std::endl;
    to_console::print_distribution_values(distribution_values, access_function);

    // Collision
    for(const auto& node : swap_information)
    {
        node_index = std::get<0>(node);
        swap_dirs = swap_directions[std::get<1>(node)];
        inout = inout_directions[std::get<2>(node)];
        std::cout << "Currently dealing with node " << node_index << std::endl;

        current_distributions = access::get_distribution_values_of(distribution_values, node_index, access_function);
        std::cout << "Distribution values of this node are " << std::endl;
        to_console::print_vector(current_distributions, current_distributions.size());

        current_velocity = macroscopic::flow_velocity(current_distributions);
        std::cout.precision(15);
        std::cout << "Flow velocity is (" << current_velocity[0] << ", " << current_velocity[1] << ")" << std::endl; 

        current_density = macroscopic::density(current_distributions);

        std::cout << "Current density is " << current_density << std::endl;
        std::cout.precision(3);
        velocities[node_index] = current_velocity;
        densities[node_index] = current_density;

        current_distributions = collision::collide_bgk(current_distributions, current_velocity, current_density);
        access::set_distribution_values_of(current_distributions, distribution_values, node_index, access_function);

        std::cout << "Distribution values after collision: " << std::endl;

        to_console::print_vector(current_distributions, current_distributions.size());
        std::cout << std::endl;
    }
    std::cout << "Distribution values after collision: " << std::endl;
    to_console::print_distribution_values(distribution_values, access_function);

    /* Update ghost nodes */
    boundary_conditions::update_velocity_input_density_output(distribution_values, access_function);
    std::cout << "Distribution values after ghost node update: " << std::endl;
    to_console::print_distribution_values(distribution_values, access_function);
    std::cout << std::endl;

    unsigned int update_node = 0;
    for(auto x = 0; x < HORIZONTAL_NODES; x = x + HORIZONTAL_NODES - 1)
    {
        for(auto y = 1; y < VERTICAL_NODES - 1; ++y)
        {
            update_node = access::get_node_index(x,y);
            current_distributions = 
                access::get_distribution_values_of(distribution_values, update_node, access_function);
                
            velocities[update_node] = macroscopic::flow_velocity(current_distributions);
            densities[update_node] = macroscopic::density(current_distributions);
        }
    }

    sim_data_tuple result{velocities, densities};

    return result;
}

/**
 * @brief Performs the combined streaming and collision step for all fluid nodes within the simulation domain.
 *        The border conditions are enforced through ghost nodes.
 */
sim_data_tuple swap_sequential::perform_swap_stream_and_collide
(
    const swap_information &swap_information,
    std::vector<double> &distribution_values,    
    const access_function access_function
)
{
    std::vector<velocity> velocities(TOTAL_NODE_COUNT, velocity{0,0});
    velocity current_velocity = {0,0};

    std::vector<double> densities(TOTAL_NODE_COUNT, -1);
    double current_density = -1;

    unsigned int node_index = 0;
    std::set<unsigned int> swap_dirs = {0,0,0};
    std::set<unsigned int> inout = {0,0,0};
    std::vector<double> current_distributions(DIRECTION_COUNT, 0);

    
    for(const auto& node : swap_information)
    {
        node_index = std::get<0>(node);
        swap_dirs = swap_directions[std::get<1>(node)];
        inout = inout_directions[std::get<2>(node)];

        // Swapping step
        swap_sequential::perform_swap_step(distribution_values, node_index, access_function, swap_dirs);

        // Restore precious order in here
        swap_sequential::restore_order(distribution_values, node_index, access_function);

        // Perform inflow and outflow, if necessary
        boundary_conditions::single_node_inout(distribution_values, node_index, inout, access_function);

        // Collision
        current_distributions = access::get_distribution_values_of(distribution_values, node_index, access_function);
        current_velocity = macroscopic::flow_velocity(current_distributions);
        current_density = macroscopic::density(current_distributions);
        velocities[node_index] = current_velocity;
        densities[node_index] = current_density;
        current_distributions = collision::collide_bgk(current_distributions, current_velocity, current_density);
        access::set_distribution_values_of(current_distributions, distribution_values, node_index, access_function);
    }

    /* Update ghost nodes */
    boundary_conditions::update_velocity_input_density_output(distribution_values, access_function);

    unsigned int update_node = 0;
    for(auto x = 0; x < HORIZONTAL_NODES; x = x + HORIZONTAL_NODES - 1)
    {
        for(auto y = 1; y < VERTICAL_NODES - 1; ++y)
        {
            update_node = access::get_node_index(x,y);
            current_distributions = 
                access::get_distribution_values_of(distribution_values, update_node, access_function);
                
            velocities[update_node] = macroscopic::flow_velocity(current_distributions);
            densities[update_node] = macroscopic::density(current_distributions);
        }
    }

    sim_data_tuple result{velocities, densities};

    return result;
}

/**
 * @brief Performs the sequential swap algorithm for the specified number of iterations.
 * 
 * @param fluid_nodes A vector containing the indices of all fluid nodes in the domain
 * @param values the vector containing the distribution values of all nodes
 * @param access_function the access function according to which the values are to be accessed
 * @param iterations this many iterations will be performed
 */
void swap_sequential::run
(  
    std::vector<unsigned int> &fluid_nodes,       
    std::vector<double> &values, 
    border_swap_information &bsi,
    access_function access_function,
    unsigned int iterations
)
{
    std::cout << "------------------------------------------------------------------------------------------------------------------------" << std::endl;
    std::cout << "Now running sequential swap algorithm for " << iterations << " iterations." << std::endl;
    std::cout << std::endl;

    std::vector<sim_data_tuple>result(
        iterations, 
        std::make_tuple(std::vector<velocity>(TOTAL_NODE_COUNT, {0,0}), std::vector<double>(TOTAL_NODE_COUNT, 0)));

    swap_information si = swap_sequential::setup_swap_information(bsi, fluid_nodes);

    for(auto time = 0; time < iterations; ++time)
    {
        std::cout << "\033[33mIteration " << time << ":\033[0m" << std::endl;
        result[time] = swap_sequential::perform_swap_stream_and_collide
        (
            si,
            values, 
            access_function
        );
        std::cout << "Finished iteration " << time << std::endl;
    }

    std::cout << std::endl;
    std::cout << std::endl;

    std::cout << "Velocity values: " << std::endl;
    std::cout << std::endl;
    for(auto i = 0; i < iterations; ++i)
    {
        std::cout << "t = " << i << std::endl;
        std::cout << "-------------------------------------------------------------------------------- " << std::endl;
        to_console::print_velocity_vector(std::get<0>(result[i]));
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
        to_console::print_vector(std::get<1>(result[i]));
        std::cout << std::endl;
    }
    std::cout << "All done, exiting simulation. " << std::endl;
   
}