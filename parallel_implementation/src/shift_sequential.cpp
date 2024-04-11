#include "../include/shift_sequential.hpp"
#include "../include/access.hpp"
#include "../include/defines.hpp"
#include "../include/boundaries.hpp"
#include "../include/collision.hpp"
#include "../include/macroscopic.hpp"
#include "../include/utils.hpp"
#include <set>
#include <iostream>
#include <stdexcept>

/**
 * @brief Retrieves an improved version of the border swap information data structure.
 *        This method does not consider inlet and outlet ghost nodes when performing bounce-back
 *        as the inserted values will be overwritten by inflow and outflow values anyways.
 * 
 * @param fluid_nodes a vector containing the indices of all fluid nodes within the simulation domain
 * @param phase_information a vector containing the phase information for every vector (true means solid)
 * @return border_swap_information see documentation of border_swap_information
 */
border_swap_information shift_sequential::retrieve_fast_border_swap_info
(
    const std::vector<unsigned int> &fluid_nodes, 
    const std::vector<bool> &phase_information
)
{
    border_swap_information result;
    for(const auto node : fluid_nodes)
    {
        std::vector<unsigned int> current_adjacencies{node};
        for(const auto direction : STREAMING_DIRECTIONS)
        {
            unsigned int current_neighbor = lbm_access::get_neighbor(node, direction);
            if(is_non_inout_ghost_node(current_neighbor, phase_information))
            {
                current_adjacencies.push_back(direction);
            }
        }
        if(current_adjacencies.size() > 1) result.push_back(current_adjacencies);
    }
    return result;
}

/**
 * @brief Performs an outstream step for all border nodes in the directions where they border non-inout ghost nodes.
 *        The distribution values will be stored in the ghost nodes in inverted order such that
 *        after this method is executed, the border nodes can be treated like regular nodes when performing an instream.
 * 
 * @param bsi 
 * @param distribution_values 
 * @param access_function 
 */
void shift_sequential::emplace_bounce_back_values
(
    const border_swap_information &bsi,
    std::vector<double> &distribution_values,
    const access_function access_function,
    const unsigned int read_offset
)
{
    std::vector<unsigned int> directions = STREAMING_DIRECTIONS;
    unsigned int current_node = 0;

    for(const auto& border_node : bsi)
    {
        current_node = border_node[0];
        directions = {border_node.begin()+1, border_node.end()};
        for(const auto direction : directions) 
        {
            distribution_values[
                access_function(lbm_access::get_neighbor(current_node + read_offset, direction), invert_direction(direction))] = 
                  distribution_values[access_function(current_node + read_offset, direction)];
        }
    }
}

/**
 * @brief Performs a combined collision and streaming step for the specified fluid node.
 * 
 * @param values the vector containing the distribution values of all nodes
 * @param access_function the access function according to which the values are to be accessed
 * @param node_index index of the fluid node whose values are to be manipulated
 * @param read_offset read offset of the current iteration (even: 0 / odd: N_e)
 * @param write_offset write offset of the current iteration (even: N_e / odd: 0)
 */


sim_data_tuple shift_sequential::perform_shift_stream_and_collide
(
    std::vector<double> &distribution_values, 
    const std::vector<unsigned int> &fluid_nodes,
    const border_swap_information &bsi,
    const access_function access_function,
    const unsigned int iteration
)
{
    unsigned int current_node = 0;
    unsigned int read_offset = 0;
    unsigned int write_offset = 0;
    std::vector<double> current_distributions(DIRECTION_COUNT, 0);
    
    std::vector<velocity> velocities(TOTAL_NODE_COUNT, velocity{0,0});
    std::vector<double> densities(TOTAL_NODE_COUNT, -1);
    std::vector<double> debug_distributions(TOTAL_NODE_COUNT, 0);

    if((iteration % 2) == 0) // Even time step, dis correct
    {
        read_offset = 0;
        write_offset = SHIFT_OFFSET;

        // Emplace bounce-back values
        for(const auto& border_node : bsi)
        {
            shift_sequential::emplace_bounce_back_values(bsi, distribution_values, access_function, read_offset);
        }

        
        for(int i = fluid_nodes.size() - 1; i >= 0; i--)
        {
            // Streaming
            current_node = fluid_nodes[i];
            shift_sequential::shift_stream(distribution_values, access_function, current_node, read_offset, write_offset);

            // Collision
            current_node = fluid_nodes[i];
            current_distributions = lbm_access::get_distribution_values_of(distribution_values, current_node + write_offset, access_function);
            velocities[current_node] = macroscopic::flow_velocity(current_distributions);
            densities[current_node] = macroscopic::density(current_distributions);
            current_distributions = collision::collide_bgk(current_distributions, velocities[current_node], densities[current_node]);
            lbm_access::set_distribution_values_of(current_distributions, distribution_values, current_node + write_offset, access_function);
        }
    }
    else
    {
        read_offset = SHIFT_OFFSET;
        write_offset = 0;

        // Emplace bounce-back values
        for(const auto& border_node : bsi)
        {
            shift_sequential::emplace_bounce_back_values(bsi, distribution_values, access_function, read_offset);
        }

        for(int i = 0; i < fluid_nodes.size(); ++i)
        {
            // Streaming
            current_node = fluid_nodes[i];
            shift_sequential::shift_stream(distribution_values, access_function, current_node, read_offset, write_offset);

            // Collision
            current_node = fluid_nodes[i];
            current_distributions = lbm_access::get_distribution_values_of(distribution_values, current_node + write_offset, access_function);
            velocities[current_node] = macroscopic::flow_velocity(current_distributions);
            densities[current_node] = macroscopic::density(current_distributions);
            current_distributions = collision::collide_bgk(current_distributions, velocities[current_node], densities[current_node]);
            lbm_access::set_distribution_values_of(current_distributions, distribution_values, current_node + write_offset, access_function);
        }
    }

    /* Update ghost nodes OFFSET ISSUE?*/
    shift_sequential::update_velocity_input_density_output(distribution_values, access_function, write_offset);
    
    unsigned int update_node = 0;
    for(auto x = 0; x < HORIZONTAL_NODES; x = x + HORIZONTAL_NODES - 1)
    {
        for(auto y = 1; y < VERTICAL_NODES - 1; ++y)
        {
            update_node = lbm_access::get_node_index(x,y);
            current_distributions = 
                lbm_access::get_distribution_values_of(distribution_values, update_node + write_offset, access_function);
                
            velocities[update_node] = macroscopic::flow_velocity(current_distributions);
            densities[update_node] = macroscopic::density(current_distributions);
        }
    }

    sim_data_tuple result{velocities, densities};
    return result;
}

/**
 * @brief Performs a combined collision and streaming step for the specified fluid node.
 * 
 * @param values the vector containing the distribution values of all nodes
 * @param access_function the access function according to which the values are to be accessed
 * @param node_index index of the fluid node whose values are to be manipulated
 * @param read_offset read offset of the current iteration (even: 0 / odd: N_e)
 * @param write_offset write offset of the current iteration (even: N_e / odd: 0)
 */


sim_data_tuple shift_sequential::perform_shift_stream_and_collide_debug
(
    std::vector<double> &distribution_values, 
    const std::vector<unsigned int> &fluid_nodes,
    const border_swap_information &bsi,
    const access_function access_function,
    const unsigned int iteration
)
{
    unsigned int current_node = 0;
    unsigned int read_offset = 0;
    unsigned int write_offset = 0;
    std::vector<double> current_distributions(DIRECTION_COUNT, 0);
    
    std::vector<velocity> velocities(TOTAL_NODE_COUNT, velocity{0,0});
    std::vector<double> densities(TOTAL_NODE_COUNT, -1);
    std::vector<double> debug_distributions(TOTAL_NODE_COUNT, 0);

    if((iteration % 2) == 0) // Even time step, dis correct
    {
        read_offset = 0;
        write_offset = SHIFT_OFFSET;

        debug_distributions = {distribution_values.begin() + (read_offset * DIRECTION_COUNT), distribution_values.end() - ((SHIFT_OFFSET - read_offset) * DIRECTION_COUNT)};
        std::cout << "Starting offset: " << read_offset * DIRECTION_COUNT << std::endl;
        std::cout << "Ending offset: " << ((SHIFT_OFFSET - read_offset) * DIRECTION_COUNT) << std::endl;
        std::cout << "Distribution values before stream and collide: " << std::endl;
        to_console::print_distribution_values(debug_distributions, access_function);
        std::cout << std::endl;

        // Emplace bounce-back values
        for(const auto& border_node : bsi)
        {
            shift_sequential::emplace_bounce_back_values(bsi, distribution_values, access_function, read_offset);
        }

        std::cout << "Distribution values after emplace bounce-back: " << std::endl;
        debug_distributions = {distribution_values.begin() + (read_offset * DIRECTION_COUNT), distribution_values.end() - ((SHIFT_OFFSET - read_offset) * DIRECTION_COUNT)};
        to_console::print_distribution_values(debug_distributions, access_function);

        std::cout << "Fluid nodes are " << std::endl;
        to_console::print_vector(fluid_nodes, HORIZONTAL_NODES - 2);

        std::cout << "Size: " << fluid_nodes.size() << std::endl;

        std::cout << "Backward iteration: " << std::endl;
        for(int i = fluid_nodes.size() - 1; i >= 0; i--)
        {
            if(i > 1000) throw std::invalid_argument("Oh no");
            std::cout << "Index " << i << ": " <<  fluid_nodes[i] << std::endl;
        }
        std::cout << std::endl;

        std::cout << "Forward iteration: " << std::endl;
        for(int i = 0; i < fluid_nodes.size(); ++i)
        {
            std::cout << "Index " << i << ": " <<  fluid_nodes[i] << std::endl;
        }
        std::cout << std::endl;

        // Streaming
        for(int i = fluid_nodes.size() - 1; i >= 0; i--)
        {
            current_node = fluid_nodes[i];
            shift_sequential::shift_stream(distribution_values, access_function, current_node, read_offset, write_offset);
        }

        std::cout << "Distribution values after streaming (properly shifted): " << std::endl;
        debug_distributions = {distribution_values.begin() + write_offset * DIRECTION_COUNT, distribution_values.end() - (SHIFT_OFFSET - write_offset) * DIRECTION_COUNT};
        to_console::print_distribution_values(debug_distributions, access_function);

        // Collision
        for(int i = fluid_nodes.size() - 1; i >= 0; --i)
        {
            current_node = fluid_nodes[i];
            std::cout << "Performing collision for node " << current_node << std::endl;
            
            std::cout << "Got distribution values " << std::endl;

            current_distributions = lbm_access::get_distribution_values_of(distribution_values, current_node + write_offset, access_function);

            to_console::print_vector(current_distributions, current_distributions.size());

            velocities[current_node] = macroscopic::flow_velocity(current_distributions);

            std::cout << "Got velocity (" << velocities[current_node][0] << ", " << velocities[current_node][1] << ")" << std::endl;

            densities[current_node] = macroscopic::density(current_distributions);

            std::cout << "Got density " << densities[current_node] << std::endl;

            std::cout << "Accessing collide_bgk " << std::endl;
            current_distributions = collision::collide_bgk(current_distributions, velocities[current_node], densities[current_node]);

            std::cout << "Setting distribution values " << std::endl;
            lbm_access::set_distribution_values_of(current_distributions, distribution_values, current_node + write_offset, access_function);
        }

        std::cout << "Distribution values after collision: " << std::endl;
        debug_distributions = {distribution_values.begin() + write_offset * DIRECTION_COUNT, distribution_values.end() - (SHIFT_OFFSET - write_offset) * DIRECTION_COUNT};
        to_console::print_distribution_values(debug_distributions, access_function);
    }
    else
    {
        read_offset = SHIFT_OFFSET;
        write_offset = 0;

        debug_distributions = {distribution_values.begin() + (read_offset * DIRECTION_COUNT), distribution_values.end() - ((SHIFT_OFFSET - read_offset) * DIRECTION_COUNT)};
        std::cout << "Starting offset: " << read_offset * DIRECTION_COUNT << std::endl;
        std::cout << "Ending offset: " << ((SHIFT_OFFSET - read_offset) * DIRECTION_COUNT) << std::endl;
        std::cout << "Distribution values before stream and collide: " << std::endl;
        to_console::print_distribution_values(debug_distributions, access_function);
        std::cout << std::endl;

        // Emplace bounce-back values
        for(const auto& border_node : bsi)
        {
            shift_sequential::emplace_bounce_back_values(bsi, distribution_values, access_function, read_offset);
        }

        std::cout << "Distribution values after emplace bounce-back: " << std::endl;
        debug_distributions = {distribution_values.begin() + (read_offset * DIRECTION_COUNT), distribution_values.end() - ((SHIFT_OFFSET - read_offset) * DIRECTION_COUNT)};
        to_console::print_distribution_values(debug_distributions, access_function);

        std::cout << "Fluid nodes are " << std::endl;
        to_console::print_vector(fluid_nodes, HORIZONTAL_NODES - 2);

        std::cout << "Size: " << fluid_nodes.size() << std::endl;

        std::cout << "Backward iteration: " << std::endl;
        for(int i = fluid_nodes.size() - 1; i >= 0; i--)
        {
            if(i > 1000) throw std::invalid_argument("Oh no");
            std::cout << "Index " << i << ": " <<  fluid_nodes[i] << std::endl;
        }
        std::cout << std::endl;

        std::cout << "Forward iteration: " << std::endl;
        for(int i = 0; i < fluid_nodes.size(); ++i)
        {
            std::cout << "Index " << i << ": " <<  fluid_nodes[i] << std::endl;
        }
        std::cout << std::endl;

        // Streaming
        for(int i = 0; i < fluid_nodes.size(); ++i)
        {
            current_node = fluid_nodes[i];
            shift_sequential::shift_stream(distribution_values, access_function, current_node, read_offset, write_offset);
        }

        std::cout << "Distribution values after streaming (properly shifted): " << std::endl;
        debug_distributions = {distribution_values.begin() + write_offset * DIRECTION_COUNT, distribution_values.end() - (SHIFT_OFFSET - write_offset) * DIRECTION_COUNT};
        to_console::print_distribution_values(debug_distributions, access_function);

        // Collision
        for(int i = 0; i < fluid_nodes.size(); ++i)
        {
            current_node = fluid_nodes[i];
            std::cout << "Performing collision for node " << current_node << std::endl;
            
            std::cout << "Got distribution values " << std::endl;

            current_distributions = lbm_access::get_distribution_values_of(distribution_values, current_node + write_offset, access_function);

            to_console::print_vector(current_distributions, current_distributions.size());

            velocities[current_node] = macroscopic::flow_velocity(current_distributions);

            std::cout << "Got velocity (" << velocities[current_node][0] << ", " << velocities[current_node][1] << ")" << std::endl;

            densities[current_node] = macroscopic::density(current_distributions);

            std::cout << "Got density " << densities[current_node] << std::endl;

            std::cout << "Accessing collide_bgk " << std::endl;
            current_distributions = collision::collide_bgk(current_distributions, velocities[current_node], densities[current_node]);

            std::cout << "Setting distribution values " << std::endl;
            lbm_access::set_distribution_values_of(current_distributions, distribution_values, current_node + write_offset, access_function);
        }

        std::cout << "Distribution values after collision: " << std::endl;
        debug_distributions = {distribution_values.begin() + write_offset * DIRECTION_COUNT, distribution_values.end() - (SHIFT_OFFSET - write_offset) * DIRECTION_COUNT};
        to_console::print_distribution_values(debug_distributions, access_function);
    }

    /* Update ghost nodes OFFSET ISSUE?*/
    shift_sequential::update_velocity_input_density_output(distribution_values, access_function, write_offset);
    
    std::cout << "Distribution values after inout update: " << std::endl;
    debug_distributions = {distribution_values.begin() + write_offset * DIRECTION_COUNT, distribution_values.end() - (SHIFT_OFFSET - write_offset) * DIRECTION_COUNT};
    to_console::print_distribution_values(debug_distributions, access_function);

    unsigned int update_node = 0;
    for(auto x = 0; x < HORIZONTAL_NODES; x = x + HORIZONTAL_NODES - 1)
    {
        for(auto y = 1; y < VERTICAL_NODES - 1; ++y)
        {
            update_node = lbm_access::get_node_index(x,y);
            current_distributions = 
                lbm_access::get_distribution_values_of(distribution_values, update_node + write_offset, access_function);
                
            velocities[update_node] = macroscopic::flow_velocity(current_distributions);
            densities[update_node] = macroscopic::density(current_distributions);
        }
    }

    sim_data_tuple result{velocities, densities};
    return result;
}

/**
 * @brief Updates the ghost nodes that represent inlet and outlet edges.
 *        When updating, a velocity border condition will be considered for the input
 *        and a density border condition for the output.
 *        The corresponding values are constants defined in "../include/"defines.hpp".
 * 
 * @param distribution_values a vector containing the distribution values of all nodes
 * @param access_function the access function used to access the distribution values
 */
void shift_sequential::update_velocity_input_density_output
(
    std::vector<double> &distribution_values, 
    const access_function access_function,
    const unsigned int offset
)
{
    // std::cout << "\033[34mEntering shift_sequential::update_velocity_input_density_output\033[0m" << std::endl;

    std::vector<double> current_dist_vals(DIRECTION_COUNT, 0);
    int current_border_node = 0;
    velocity v = INLET_VELOCITY;
    double density = 0;

    // // Correction of bordering velocity
    // current_border_node = lbm_access::get_node_index(HORIZONTAL_NODES - 2,0) + offset;
    // current_dist_vals = lbm_access::get_distribution_values_of(distribution_values, lbm_access::get_neighbor(current_border_node, 7), access_function);
    // lbm_access::set_distribution_values_of(current_dist_vals, distribution_values, current_border_node, access_function);

    // current_border_node = lbm_access::get_node_index(HORIZONTAL_NODES - 2,VERTICAL_NODES - 1) + offset;
    // current_dist_vals = lbm_access::get_distribution_values_of(distribution_values, lbm_access::get_neighbor(current_border_node, 1), access_function);
    // lbm_access::set_distribution_values_of(current_dist_vals, distribution_values, current_border_node, access_function);


    for(auto y = 0; y < VERTICAL_NODES; ++y)
    {
        // Update inlets
        // std::cout << "Currently at node with coords (0, " << y << "); index is " << access::get_node_index(0,y) << " (with offset added: " << access::get_node_index(0,y) + offset << ")" << std::endl;
        current_border_node = lbm_access::get_node_index(0,y) + offset;
        v = INLET_VELOCITY;
        density = INLET_DENSITY;
        current_dist_vals = maxwell_boltzmann_distribution(v, density);
        lbm_access::set_distribution_values_of
        (
            current_dist_vals,
            distribution_values,
            current_border_node,
            access_function
        );

        // std::cout << std::endl;

        // Update outlets
        // std::cout << "Currently at node with coords (" << HORIZONTAL_NODES - 1 << ", " << y << "); index is " 
        //          << access::get_node_index(HORIZONTAL_NODES - 1,y) << " (with offset added: " 
        //          << access::get_node_index(HORIZONTAL_NODES - 1,y) + offset << ")" << std::endl;

        current_border_node = lbm_access::get_node_index(HORIZONTAL_NODES - 1,y) + offset;

        v = macroscopic::flow_velocity(lbm_access::get_distribution_values_of(distribution_values, lbm_access::get_neighbor(current_border_node, 3), access_function));

        // std::cout << "\tVelocity is (" << v[0] << ", " << v[1] << "), got distribution values from neighbor node with index " << access::get_neighbor(current_border_node, 3) << std::endl;
        // std::cout << "This neighbor has the distribution values " << std::endl;
        // to_console::print_vector(access::get_distribution_values_of(distribution_values, access::get_neighbor(current_border_node, 3), access_function), 100);

        density = OUTLET_DENSITY;

        // std::cout << "\tDensity is " << density << std::endl;

        current_dist_vals = maxwell_boltzmann_distribution(v, density);
        // std::cout << "\tCalculated Maxwell Boltzmann distribution " << std::endl;
        // to_console::print_vector(current_dist_vals, current_dist_vals.size());
        lbm_access::set_distribution_values_of
        (
            current_dist_vals,
            distribution_values,
            current_border_node,
            access_function
        );
        // std::cout << std::endl;
    }
}

/**
 * @brief Performs the sequential shift algorithm for the specified number of iterations.
 * 
 * @param fluid_nodes A vector containing the indices of all fluid nodes in the domain
 * @param values the vector containing the distribution values of all nodes
 * @param access_function the access function according to which the values are to be accessed
 * @param iterations this many iterations will be performed
 */
void shift_sequential::run
(  
    std::vector<unsigned int> &fluid_nodes,       
    std::vector<double> &values, 
    border_swap_information &bsi,
    access_function access_function,
    unsigned int iterations
)
{
    std::cout << "------------------------------------------------------------------------------------------------------------------------" << std::endl;
    std::cout << "Now running sequential shift algorithm for " << iterations << " iterations." << std::endl;
    std::cout << std::endl;

    std::vector<sim_data_tuple>result(
        iterations, 
        std::make_tuple(std::vector<velocity>(TOTAL_NODE_COUNT, {0,0}), std::vector<double>(TOTAL_NODE_COUNT, 0)));

    for(auto time = 0; time < iterations; ++time)
    {
        std::cout << "\033[33mIteration " << time << ":\033[0m" << std::endl;
        result[time] = shift_sequential::perform_shift_stream_and_collide(values, fluid_nodes, bsi, access_function, time);
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

/**
 * @brief Create an example domain for testing purposes. The domain is a rectangle with
 *        dimensions specified in the defines file where the outermost nodes are ghost nodes.
 *        The upper and lower ghost nodes are solid whereas the leftmost and rightmost columns are fluid
 *        nodes that mark the inlet and outlet respectively.
 *        Notice that all data will be written to the parameters which are assumed to be empty initially.
 * 
 * @param distribution_values a vector containing all distribution values.
 * @param nodes a vector containing all node indices, including those of solid nodes and ghost nodes.
 * @param fluid_nodes a vector containing the indices of all fluid nodes.
 * @param phase_information a vector containing the phase information of all nodes where true means solid.
 * @param access_function the domain will be prepared for access with this access function.
 */
void shift_sequential::setup_example_domain
(
    std::vector<double> &distribution_values,
    std::vector<unsigned int> &nodes,
    std::vector<unsigned int> &fluid_nodes,
    std::vector<bool> &phase_information,
    border_swap_information &swap_info,
    access_function access_function
)
{
    std::cout << "Setting up example domain." << std::endl;

    /* Set up distribution values */

    // std::vector<double> values_0 = {0.00,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08};
    // std::vector<double> values_1 = {-0.00,0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008}; 

    std::vector<double> values = maxwell_boltzmann_distribution(VELOCITY_VECTORS.at(4), 1);

    for(auto i = 0; i < TOTAL_NODE_COUNT; ++i)
    {
        // if(i % 2) distribution_values.insert(distribution_values.end(), values_0.begin(), values_0.end());
        // else distribution_values.insert(distribution_values.end(), values_1.begin(), values_1.end());

        distribution_values.insert(distribution_values.end(), values.begin(), values.end());
    }

    for(auto i = 0; i < SHIFT_OFFSET; ++i)
    {
        distribution_values.insert(distribution_values.end(), {0,0,0,0,0,0,0,0,0});
    }
    
    boundary_conditions::initialize_inout(distribution_values, access_function);

    std::cout << "All distribution values were set." << std::endl;

    /* Set all nodes for direct access */
    for(auto i = 0; i < TOTAL_NODE_COUNT; ++i)
    {
        nodes.push_back(i);
    }

    /* Set up vector containing fluid nodes within the simulation domain. */
    for(auto it = nodes.begin() + HORIZONTAL_NODES; it < nodes.end() - HORIZONTAL_NODES; ++it)
    {
        if(((*it % HORIZONTAL_NODES) != 0) && ((*it % HORIZONTAL_NODES) != (HORIZONTAL_NODES - 1))) fluid_nodes.push_back(*it);
    }
    
    /* Phase information vector */
    phase_information.assign(TOTAL_NODE_COUNT, false);
    for(auto x = 0; x < HORIZONTAL_NODES; ++x)
    {
        phase_information[lbm_access::get_node_index(x,0)] = true;
        phase_information[lbm_access::get_node_index(x,VERTICAL_NODES - 1)] = true;
    }

    /* Set up border swap information */
    swap_info = shift_sequential::retrieve_fast_border_swap_info(fluid_nodes, phase_information);
}