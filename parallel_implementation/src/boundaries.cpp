#include "../include/boundaries.hpp"
#include "../include/macroscopic.hpp"
#include "../include/access.hpp"
#include "../include/utils.hpp"
#include <iostream>
#include <hpx/format.hpp>
#include <hpx/future.hpp>
#include <hpx/algorithm.hpp>
#include <hpx/execution.hpp>
#include <hpx/iostream.hpp>

/**
 * @brief Returns whether the node with the specified index is located at the edge of the simulation domain.
 *        This is the case for any of the following coordinates:
 *        - (1, y)
 *        - (HORIZONTAL_NODES - 2, y)
 *        - (x, 1)
 *        - (x, VERTICAL_NODES - 2)
 *        with suitable x and y.
 * 
 * @param node_index the index of the node in question
 * @return whether or not the node is an edge node
 *         
 */
bool is_edge_node(unsigned int node_index)
{
    std::tuple<unsigned int, unsigned int> coordinates = lbm_access::get_node_coordinates(node_index);
    unsigned int x = std::get<0>(coordinates);
    unsigned int y = std::get<1>(coordinates);
    bool result = ((x == 1) || (x == (HORIZONTAL_NODES - 2))) && ((y == 1) || (y == (VERTICAL_NODES - 2)));
    return result;
}

/**
 * @brief Returns whether the node with the specified index is a ghost node.
 *        This is the case for any of the following:
 *        
 *        With suitable x and y, the node coordinates are any of
 *        - (0, y)
 *        - (HORIZONTAL_NODES - 1, y)
 *        - (x, 0)
 *        - (x, VERTICAL_NODES - 1)
 *        
 *        or
 * 
 *        The node with the specified index is solid.
 * 
 * @param node_index the index of the node in question
 * @param phase_information a vector containing the phase information for all nodes of the lattice
 */
bool is_ghost_node(unsigned int node_index, const std::vector<bool> &phase_information)
{
    std::tuple<unsigned int, unsigned int> coordinates = lbm_access::get_node_coordinates(node_index);
    unsigned int x = std::get<0>(coordinates);
    unsigned int y = std::get<1>(coordinates);
    bool is_outer_node = ((x == 0) || (x == (HORIZONTAL_NODES - 1))) || ((y == 0) || (y == (VERTICAL_NODES - 1)));
    return is_outer_node  | phase_information[node_index];
}

/**
 * @brief Retrieves the border swap information for all fluid nodes within the simulation domain based on the 
 *        phase information of all nodes. Notice that all fluid nodes on the edges of the simulation domain will 
 *        automatically become border nodes.
 * 
 * @param fluid_nodes a vector containing the indices of all fluid nodes within the simulation domain
 * @param phase_information a vector containing the phase information of ALL nodes 
 * @return see documentation of border_swap_information
 */
border_swap_information bounce_back::retrieve_border_swap_information
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
            if(is_ghost_node(current_neighbor, phase_information))
            {
                current_adjacencies.push_back(direction);
            }
        }
        if(current_adjacencies.size() > 1) result.push_back(current_adjacencies);
    }
    return result;
}

/**
 * @brief Retrieves an improved version of the border swap information data structure.
 *        This method does not consider inlet and outlet ghost nodes when performing bounce-back
 *        as the inserted values will be overwritten by inflow and outflow values anyways.
 * 
 * @param fluid_nodes a vector containing the indices of all fluid nodes within the simulation domain
 * @param phase_information a vector containing the phase information for every vector (true means solid)
 * @return border_swap_information see documentation of border_swap_information
 */
border_swap_information bounce_back::retrieve_fast_border_swap_info
(
    const std::vector<unsigned int> &fluid_nodes, 
    const std::vector<bool> &phase_information
)
{
    std::vector<unsigned int> current_adjacencies;
    border_swap_information result;
    for(const auto node : fluid_nodes)
    {
        current_adjacencies = {node};
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
 * @param bsi a border_swap_information generated by retrieve_fast_border_swap_info
 * @param distribution_values a vector containing the distribution values of all nodes
 * @param access_function the access function used to access the distribution values
 * @param read_offset offset for shift algorithm, leave zero for all other algorithms
 */
void bounce_back::emplace_bounce_back_values
(
    const border_swap_information &bsi,
    std::vector<double> &distribution_values,
    const access_function access_function,
    const unsigned int read_offset
)
{
    for(auto bsi_iterator = bsi.begin(); bsi_iterator < bsi.end(); ++bsi_iterator)
    {
        for(auto direction_iterator = (*bsi_iterator).begin()+1; direction_iterator < (*bsi_iterator).end(); ++direction_iterator) 
        {
            distribution_values[
                access_function(lbm_access::get_neighbor((*bsi_iterator)[0] + read_offset, *direction_iterator), invert_direction(*direction_iterator))] = 
                  distribution_values[access_function((*bsi_iterator)[0] + read_offset, *direction_iterator)];
        }
    }
}

/**
 * @brief Performs an outstream step for all border nodes in the directions where they border non-inout ghost nodes.
 *        The distribution values will be stored in the ghost nodes in inverted order such that
 *        after this method is executed, the border nodes can be treated like regular nodes when performing an instream.
 * 
 * @param bsi a border_swap_information generated by retrieve_fast_border_swap_info
 * @param distribution_values a vector containing the distribution values of all nodes
 * @param access_function the access function used to access the distribution values
 * @param read_offset offset for shift algorithm, leave zero for all other algorithms
 */
void bounce_back::emplace_bounce_back_values_parallel
(
    const border_swap_information &bsi,
    std::vector<double> &distribution_values,
    const access_function access_function,
    const unsigned int read_offset
)
{
    hpx::for_each
    (
        hpx::execution::par, 
        bsi.begin(), 
        bsi.end(), 
        [&](const std::vector<unsigned int>& fluid_node)
        {
            for(auto direction_iterator = fluid_node.begin()+1; direction_iterator < fluid_node.end(); ++direction_iterator) 
            {
                distribution_values[
                    access_function(lbm_access::get_neighbor(fluid_node[0] + read_offset, *direction_iterator), invert_direction(*direction_iterator))] = 
                    distribution_values[access_function(fluid_node[0] + read_offset, *direction_iterator)];
            }
        }
    );
}

/**
 * @brief Determines the directions in which the value will be reflected from the opposite direction.
 * 
 * @param current_border_info an entry of a border_swap_information object
 * @return a set containing all bounce back directions
 */
std::set<unsigned int> bounce_back::determine_bounce_back_directions
(
    const std::vector<unsigned int> &current_border_info
)
{
    std::vector<unsigned int> inverted;
    for(auto j = current_border_info.begin() + 1; j < current_border_info.end(); ++j)
    {
        inverted.push_back(invert_direction(*j));
    }
    std::set<unsigned int> remaining_dirs = {inverted.begin(), inverted.end()};
    return remaining_dirs;  
}

/**
 * @brief Performs a halfway bounce-back streaming update for all fluid nodes within the simulation domain.
 *        This version utilizes the ghost nodes bordering a boundary node. It is intended for use with
 *        the two-step, swap and shift algorithms.
 * 
 * @param ba see documentation of border_adjacency
 * @param distribution_values a vector containing the distribution values of all nodes
 * @param access_function the access function used to access the distribution values
 */
void bounce_back::perform_boundary_update
(
    const border_swap_information &bsi,
    std::vector<double> &distribution_values, 
    const access_function access_function
)
{
    std::vector<double> current_dist_vals;
    std::vector<std::tuple<unsigned int, unsigned int>> swap_partners;
    std::set<unsigned int> remaining_dirs{STREAMING_DIRECTIONS.begin(), STREAMING_DIRECTIONS.end()};
    int current_border_node = 0;
    for(const auto& current : bsi)
    {
        //std::cout << "Received bsi " << std::endl;
        //to_console::print_vector(current, current.size());
        current_border_node = current[0];
        remaining_dirs = bounce_back::determine_bounce_back_directions(current);
        //std::cout << "Remaining dirs: " << std::endl;
        //to_console::print_set(remaining_dirs);
        //std::cout << std::endl;

        for(const auto direction : remaining_dirs)
        {
            distribution_values[access_function(current_border_node, direction)] = 
            distribution_values[access_function(lbm_access::get_neighbor(current_border_node, invert_direction(direction)), invert_direction(direction))];
        }
    }
}

/**
 * @brief Modified version of the halfway bounce-back streaming update for all fluid nodes 
 *        within the simulation domain. Instead of using information stored in ghost nodes, 
 *        This allows for a convenient unification of the streaming and collision step for
 *        the two-lattice algorithm.
 * 
 * @param bsi see documentation of border_swap_information
 * @param source the distribution values will be read from this vector
 * @param destination the updated distribution values will be written to this vector
 * @param access_function the access function used to access the distribution values
 */
void bounce_back::perform_early_boundary_update
(
    const border_swap_information &bsi,
    const std::vector<double> &source, 
    std::vector<double> &destination, 
    const access_function access_function
)
{
    std::vector<double> current_dist_vals;
    std::set<unsigned int> remaining_dirs{STREAMING_DIRECTIONS.begin(), STREAMING_DIRECTIONS.end()};
    int current_border_node = 0;
    for(const auto& current : bsi)
    {
        current_border_node = current[0];
        remaining_dirs = bounce_back::determine_bounce_back_directions(current);
        
        for(const auto direction : remaining_dirs)
        {
            destination[access_function(current_border_node, direction)] = 
            source[access_function(current_border_node, invert_direction(direction))];
        }
    }
}

/**
 * @brief Modified version of the halfway bounce-back streaming update for all fluid nodes 
 *        within the simulation domain. Instead of using information stored in ghost nodes, 
 *        This allows for a convenient unification of the streaming and collision step for
 *        the two-lattice algorithm.
 *        This variant will update inlets and outlets according to the specified velocity and density.
 * 
 * @param bsi see documentation of border_swap_information
 * @param distribution_values a vector containing the distribution values of all nodes
 * @param access_function the access function used to access the distribution values
 */
void boundary_conditions::perform_inout_boundary_update
(
    const border_swap_information &bsi,
    std::vector<double> &distribution_values, 
    const access_function access_function
)
{
    std::vector<double> current_dist_vals(DIRECTION_COUNT, 0);
    std::set<unsigned int> remaining_dirs{STREAMING_DIRECTIONS.begin(), STREAMING_DIRECTIONS.end()};
    int current_border_node = 0;
    velocity v = INLET_VELOCITY;

    for(const auto& current : bsi)
    {
        current_border_node = current[0];
        if(std::get<0>(lbm_access::get_node_coordinates(current_border_node)) == 1) // Inlet node
        {
            v = INLET_VELOCITY;
            current_dist_vals = maxwell_boltzmann_distribution(v, INLET_DENSITY);
            lbm_access::set_distribution_values_of
            (
                current_dist_vals,
                distribution_values,
                current_border_node,
                access_function
            );
        }
        else if(std::get<0>(lbm_access::get_node_coordinates(current_border_node)) == HORIZONTAL_NODES - 1) // Outlet node
        {
            v = OUTLET_VELOCITY;
            current_dist_vals = maxwell_boltzmann_distribution(v, OUTLET_DENSITY);
            lbm_access::set_distribution_values_of
            (
                current_dist_vals,
                distribution_values,
                current_border_node,
                access_function
            );
        }
        else
        {
            remaining_dirs = bounce_back::determine_bounce_back_directions(current);
            for(auto direction : remaining_dirs)
            {
                distribution_values[access_function(current_border_node, direction)] = 
                distribution_values[access_function(current_border_node, invert_direction(direction))];
            }
        }

    }
}

/**
 * @brief Updates the ghost nodes that represent inlet and outlet edges.
 *        When updating, a velocity border condition will be considered for both the input and the output.
 *        The corresponding values are constants defined in "../include/"defines.hpp".
 * 
 * @param distribution_values a vector containing the distribution values of all nodes
 * @param velocities a vector containing the velocities of all nodes
 * @param densities a vector containing the densities of all nodes
 * @param access_function the access function used to access the distribution values
 */
void boundary_conditions::update_velocity_input_velocity_output
(
    std::vector<double> &distribution_values,
    std::vector<velocity> &velocities,
    std::vector<double> &densities, 
    const access_function access_function
)
{
    std::vector<double> current_dist_vals(DIRECTION_COUNT, 0);
    int current_border_node = 0;
    velocity in = INLET_VELOCITY;
    velocity out = OUTLET_VELOCITY;
    std::vector<velocity> inlet = velocity_profiles::ideal_laminary(in);
    std::vector<velocity> outlet = velocity_profiles::seventh_rule_turbulent(out);
    double density = 0;

    for(auto y = 1; y < VERTICAL_NODES - 1; ++y)
    {
        // Update inlets
        current_border_node = lbm_access::get_node_index(0,y);
        density = macroscopic::density(lbm_access::get_distribution_values_of(distribution_values, lbm_access::get_neighbor(current_border_node, 5), access_function));
        density = INLET_DENSITY + (INLET_DENSITY - density);
        current_dist_vals = maxwell_boltzmann_distribution(inlet[y-1], density);
        lbm_access::set_distribution_values_of
        (
            current_dist_vals,
            distribution_values,
            current_border_node,
            access_function
        );
        velocities[current_border_node] = inlet[y-1];
        densities[current_border_node] = density;

        // Update outlets
        current_border_node = lbm_access::get_node_index(HORIZONTAL_NODES - 1,y);
        density = macroscopic::density(lbm_access::get_distribution_values_of(distribution_values, lbm_access::get_neighbor(current_border_node, 3), access_function));
        density = OUTLET_DENSITY + (OUTLET_DENSITY - density);
        current_dist_vals = maxwell_boltzmann_distribution(outlet[y-1], density);
        lbm_access::set_distribution_values_of
        (
            current_dist_vals,
            distribution_values,
            current_border_node,
            access_function
        );
        velocities[current_border_node] = outlet[y-1];
        densities[current_border_node] = density;
    }
}

/**
 * @brief Updates the ghost nodes that represent inlet and outlet edges.
 *        When updating, a velocity border condition will be considered for the input
 *        and a density border condition for the output.
 *        The inlet velocity is constant throughout all inlet nodes whereas the outlet nodes
 *        all have the specified density.
 *        The corresponding values are constants defined in "../include/"defines.hpp".
 * 
 * @param distribution_values a vector containing the distribution values of all nodes
 * @param velocities a vector containing the velocities of all nodes
 * @param densities a vector containing the densities of all nodes
 * @param access_function the access function used to access the distribution values
 */
void boundary_conditions::update_velocity_input_density_output
(
    std::vector<double> &distribution_values,
    std::vector<velocity> &velocities,
    std::vector<double> &densities, 
    const access_function access_function
)
{
    std::vector<double> current_dist_vals(DIRECTION_COUNT, 0);
    int current_border_node = 0;
    velocity v = INLET_VELOCITY;
    double density = 0;

    for(auto y = 0; y < VERTICAL_NODES; ++y)
    {
        // Update inlets
        current_border_node = lbm_access::get_node_index(0,y);
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
        velocities[current_border_node] = v;
        densities[current_border_node] = density;

        // Update outlets
        current_border_node = lbm_access::get_node_index(HORIZONTAL_NODES - 1,y);
        v = macroscopic::flow_velocity(lbm_access::get_distribution_values_of(distribution_values, lbm_access::get_neighbor(current_border_node, 3), access_function));
        density = OUTLET_DENSITY;
        current_dist_vals = maxwell_boltzmann_distribution(v, density);
        lbm_access::set_distribution_values_of
        (
            current_dist_vals,
            distribution_values,
            current_border_node,
            access_function
        );
        velocities[current_border_node] = v;
        densities[current_border_node] = density;
    }
}

/**
 * @brief Updates the ghost nodes that represent inlet and outlet edges.
 *        When updating, a velocity border condition will be considered for the input
 *        and a density border condition for the output.
 *        The inlet velocity is constant throughout all inlet nodes whereas the outlet nodes
 *        all have the specified density.
 *        The corresponding values are constants defined in "../include/"defines.hpp".
 * 
 * @param distribution_values a vector containing the distribution values of all nodes
 * @param velocities a vector containing the velocities of all nodes
 * @param densities a vector containing the densities of all nodes
 * @param access_function the access function used to access the distribution values
 */
void boundary_conditions::update_velocity_input_density_output_parallel
(
    std::vector<double> &distribution_values,
    std::vector<velocity> &velocities,
    std::vector<double> &densities, 
    const access_function access_function
)
{
    hpx::experimental::for_loop(
        hpx::execution::par, 0, VERTICAL_NODES - 1,
        [&distribution_values, &velocities, &densities, access_function](int y)
        {
        // Update inlets
        int current_border_node = lbm_access::get_node_index(0,y);
        velocity v = INLET_VELOCITY;
        double density = INLET_DENSITY;
        std::vector<double> current_dist_vals = maxwell_boltzmann_distribution(v, density);
        lbm_access::set_distribution_values_of
        (
            current_dist_vals,
            distribution_values,
            current_border_node,
            access_function
        );
        velocities[current_border_node] = v;
        densities[current_border_node] = density;

        // Update outlets
        current_border_node = lbm_access::get_node_index(HORIZONTAL_NODES - 1,y);
        v = macroscopic::flow_velocity(lbm_access::get_distribution_values_of(distribution_values, lbm_access::get_neighbor(current_border_node, 3), access_function));
        density = OUTLET_DENSITY;
        current_dist_vals = maxwell_boltzmann_distribution(v, density);
        lbm_access::set_distribution_values_of
        (
            current_dist_vals,
            distribution_values,
            current_border_node,
            access_function
        );
        velocities[current_border_node] = v;
        densities[current_border_node] = density;
        });
}

/**
 * @brief Updates the ghost nodes that represent inlet and outlet edges.
 *        When updating, a density border condition will be considered for both the input and the output.
 *        The corresponding values are constants defined in "../include/"defines.hpp".
 * 
 * @param distribution_values a vector containing the distribution values of all nodes
 * @param velocities a vector containing the velocities of all nodes
 * @param densities a vector containing the densities of all nodes
 * @param access_function the access function used to access the distribution values
 */
void boundary_conditions::update_density_input_density_output
(
    std::vector<double> &distribution_values, 
    std::vector<velocity> &velocities,
    std::vector<double> &densities, 
    const access_function access_function
)
{
    std::vector<double> current_dist_vals(DIRECTION_COUNT, 0);
    int current_border_node = 0;
    velocity v = INLET_VELOCITY;
    double density = 0;

    for(auto y = 0; y < VERTICAL_NODES; ++y)
    {
        // Update inlets
        current_border_node = lbm_access::get_node_index(0,y);
        v = {0,0};
        density = INLET_DENSITY;
        current_dist_vals = maxwell_boltzmann_distribution(v, density);
        lbm_access::set_distribution_values_of
        (
            current_dist_vals,
            distribution_values,
            current_border_node,
            access_function
        );
        velocities[current_border_node] = v;
        densities[current_border_node] = density;

        // Update outlets
        current_border_node = lbm_access::get_node_index(HORIZONTAL_NODES - 1,y);
        v = macroscopic::flow_velocity(lbm_access::get_distribution_values_of(distribution_values, lbm_access::get_neighbor(current_border_node, 3), access_function));
        density = OUTLET_DENSITY;
        current_dist_vals = maxwell_boltzmann_distribution(v, density);
        lbm_access::set_distribution_values_of
        (
            current_dist_vals,
            distribution_values,
            current_border_node,
            access_function
        );
        velocities[current_border_node] = v;
        densities[current_border_node] = density;
    }
}

/**
 * @brief Initializes all inlet and outlet nodes with their corresponding initial values.
 *        The corresponding values are constants defined in "../include/"defines.hpp".
 * 
 * @param distribution_values a vector containing the distribution values of all nodes
 * @param access_function the access function used to access the distribution values
 */
void boundary_conditions::initialize_inout
(
    std::vector<double> &distribution_values, 
    const access_function access_function
)
{
    std::vector<double> current_dist_vals(DIRECTION_COUNT, 0);
    int current_border_node = 0;
    velocity v = INLET_VELOCITY;
    double density = 0;

    for(auto y = 0; y < VERTICAL_NODES; ++y)
    {
        // Update inlets
        current_border_node = lbm_access::get_node_index(0,y);
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

        // Update outlets
        current_border_node = lbm_access::get_node_index(HORIZONTAL_NODES - 1,y);
        v = OUTLET_VELOCITY;
        density = OUTLET_DENSITY;
        current_dist_vals = maxwell_boltzmann_distribution(v, density);
        lbm_access::set_distribution_values_of
        (
            current_dist_vals,
            distribution_values,
            current_border_node,
            access_function
        );
    }
}

/**
 * @brief Realizes inflow and outflow by an inward stream of each border node.
 *        This method is intended for use with two-step, swap and shift algorithms.
 * 
 * @param distribution_values a vector containing the distribution values of all nodes
 * @param access_function the access function used to access the distribution values
 */
void boundary_conditions::ghost_stream_inout
(
    std::vector<double> &distribution_values, 
    const access_function access_function
)
{
    std::set<unsigned int> inflow_instream_dirs{2,5,8};
    std::set<unsigned int> outflow_instream_dirs{0,3,6};
    int current_border_node = 0;

    for(auto y = 1; y < VERTICAL_NODES - 1; ++y)
    {
        // Update inlets
        current_border_node = lbm_access::get_node_index(1,y);
        for(const auto direction : inflow_instream_dirs)
        {
            distribution_values[access_function(current_border_node, direction)] = 
                distribution_values[access_function(lbm_access::get_neighbor(current_border_node, invert_direction(direction)), direction)];
        }

        // Update outlets
        current_border_node = lbm_access::get_node_index(HORIZONTAL_NODES - 2,y);
        for(const auto direction : outflow_instream_dirs)
        {
            distribution_values[access_function(current_border_node, direction)] = 
                distribution_values[access_function(lbm_access::get_neighbor(current_border_node, invert_direction(direction)), direction)];
        }
    }
}

/**
 * @brief Computes a laminary velocity profile for inlet or outlet nodes.
 * 
 * @param u the mean velocity of the profile
 * @return std::vector<velocity> a vector containing the velocity values for the inlet or outlet nodes.
 */
std::vector<velocity> velocity_profiles::ideal_laminary(velocity &u)
{
    std::vector<velocity> result;
    double middle_line = (VERTICAL_NODES)/2.0f;
    double radius = (VERTICAL_NODES - 2)/2.0f;
    for(auto y = 1; y < VERTICAL_NODES; ++y)
    {
        result.push_back({2 * INLET_VELOCITY[0] * (1 - pow((y + 0.5f - middle_line)/radius,2)) ,0});
    }
    return result;
}

/**
 * @brief Computes a turbulent velocity profile for inlet or outlet nodes using the rule of the seventh.
 * 
 * @param u the mean velocity of the profile
 * @return std::vector<velocity> a vector containing the velocity values for the inlet or outlet nodes.
 */
std::vector<velocity> velocity_profiles::seventh_rule_turbulent(velocity &u)
{
    std::vector<velocity> result;
    double middle_line = (VERTICAL_NODES)/2.0f;
    double radius = (VERTICAL_NODES - 2)/2.0f;
    for(auto y = 1; y < VERTICAL_NODES; ++y)
    {
        result.push_back({1.1f * OUTLET_VELOCITY[0] * (1 - pow((abs(y + 0.5f - middle_line))/radius,7)) ,0});
    }
    return result;
}