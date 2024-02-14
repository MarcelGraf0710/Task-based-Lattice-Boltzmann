#include "../include/boundaries.hpp"
#include "../include/macroscopic.hpp"
#include "../include/access.hpp"
#include "../include/utils.hpp"
#include <iostream>

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
    std::tuple<unsigned int, unsigned int> coordinates = access::get_node_coordinates(node_index);
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
bool is_ghost_node(unsigned int node_index, std::vector<bool> &phase_information)
{
    std::tuple<unsigned int, unsigned int> coordinates = access::get_node_coordinates(node_index);
    unsigned int x = std::get<0>(coordinates);
    unsigned int y = std::get<1>(coordinates);
    bool is_outer_node = ((x == 0) || (x == (HORIZONTAL_NODES - 1))) || ((y == 0) || (y == (VERTICAL_NODES - 1)));
    return is_outer_node  | phase_information[node_index];
}

/**
 * @brief Returns a vector containing all fluid non-border nodes within the simulation domain.
 * 
 * @param fluid_nodes a vector containing the indices of all fluid nodes within the simulation domain
 * @param ba see documentation of border_adjacency
 */
std::vector<unsigned int> get_non_border_nodes
(
    std::vector<unsigned int> &fluid_nodes,
    border_adjacency &ba
)
{
    std::vector<unsigned int> result;
    border_adjacency::iterator next_border_node = ba.begin();
    for(auto i = 0; i < fluid_nodes.size(); ++i)
    {
        if((fluid_nodes[i] == std::get<0>((*next_border_node)[0])) 
            & (next_border_node < ba.end() - 1))
        {
            ++next_border_node;
        }
        else
        {
            result.push_back(fluid_nodes[i]);
        }
    }
    return result;
}

/**
 * @brief Retrieves the border adjacencies for all fluid nodes within the simulation domain based on 
 *        the phase information of all nodes. Notice that all fluid nodes on the edges of the simulation 
 *        domain will automatically become border nodes.
 * 
 * @param fluid_nodes a vector containing the indices of all fluid nodes within the simulation domain
 * @param phase_information a vector containing the phase information of ALL nodes 
 * @return the border adjancency relations of each node as it is required by bounce-back boundary treatment
 *   
 */
border_adjacency bounce_back::retrieve_border_adjacencies
(        
    std::vector<unsigned int> &fluid_nodes, 
    std::vector<bool> &phase_information
)
{
    border_adjacency result;
    std::vector<std::tuple<unsigned int, unsigned int>> current_adjacencies;
    unsigned int current_neighbor;

    for(auto node : fluid_nodes)
    {
        current_adjacencies = {std::make_tuple(node, 4)};
        for(auto direction : streaming_directions)
        {
            current_neighbor = access::get_neighbor(node, direction);
            if(is_ghost_node(current_neighbor, phase_information))
            {
                current_adjacencies.push_back(std::make_tuple(current_neighbor, direction));
            }
        }
        if(current_adjacencies.size() > 1) result.push_back(current_adjacencies);
    }
    return result;
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
    std::vector<unsigned int> &fluid_nodes, 
    std::vector<bool> &phase_information
)
{
    border_swap_information result;
    for(auto node : fluid_nodes)
    {
        std::vector<unsigned int> current_adjacencies{node};
        for(auto direction : streaming_directions)
        {
            unsigned int current_neighbor = access::get_neighbor(node, direction);
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
     * @brief Determines the directions in which the value will be reflected from the opposite direction.
     * 
     * @param current_border_info an entry of a border_swap_information object
     * @return a set containing all bounce back directions
     */
std::set<unsigned int> bounce_back::determine_bounce_back_directions
(
    std::vector<unsigned int> &current_border_info
)
{
    //std::cout << "Determining bounce back directions for ";
    //to_console::print_vector(current_border_info, 10);
    
    std::vector<unsigned int> inverted;
    for(auto j = current_border_info.begin() + 1; j < current_border_info.end(); ++j)
    {
        inverted.push_back(invert_direction(*j));
    }
    //std::cout << "Inverted dirs are ";
    //to_console::print_vector(inverted, 10);
    std::set<unsigned int> remaining_dirs = {inverted.begin(), inverted.end()};
    // for (auto i = current_border_info.begin() + 1; i < current_border_info.end(); ++i)
    // {
    //     remaining_dirs.erase(*i);
    // }
    //std::cout << "Returning remaining dirs ";
    //to_console::print_set(remaining_dirs);
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
    border_adjacency &ba,
    std::vector<double> &distribution_values, 
    access_function access_function
)
{
    int current_border_node = 0;
    std::vector<std::tuple<unsigned int, unsigned int>> neighbors;
    unsigned int current_dir = 0;

    for(auto current_border_node_information : ba)
    {
        current_border_node = std::get<0>(current_border_node_information[0]);
        neighbors = 
        {
            current_border_node_information.begin() + 1, 
            current_border_node_information.end()
        };
        for(auto neighbor : neighbors)
        {
            current_dir = invert_direction(std::get<1>(neighbor));
            distribution_values[access_function(current_border_node, current_dir)] = 
            distribution_values[access_function(std::get<0>(neighbor), current_dir)];
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
 * @param distribution_values a vector containing the distribution values of all nodes
 * @param access_function the access function used to access the distribution values
 */
void bounce_back::perform_early_boundary_update
(
    border_swap_information &bsi,
    std::vector<double> &source, // NEW!!
    std::vector<double> &destination, // RENAMED
    access_function access_function
)
{
    std::vector<double> current_dist_vals;
    std::set<unsigned int> remaining_dirs{streaming_directions.begin(), streaming_directions.end()};
    int current_border_node = 0;
    for(auto current : bsi)
    {
        current_border_node = current[0];
        remaining_dirs = bounce_back::determine_bounce_back_directions(current);
        // if (current_border_node == 8 || current_border_node == 12 || current_border_node == 22 || current_border_node == 26)
        // {
        //     std::cout << "Currently treating node " << current_border_node << std::endl;
        //     std::cout << "\t with distribution values \t";
        //     current_dist_vals = access::get_distribution_values_of(destination, current_border_node, access_function);
        //     to_console::print_vector(current_dist_vals, 10);
        //     std::cout << "\t and directions \t";
        //     to_console::print_set(remaining_dirs);
        // }
        //std::cout << "Currently treating node " << current_border_node << std::endl;
        //std::cout << "\t with distribution values \t";
        //current_dist_vals = access::get_distribution_values_of(distribution_values, current_border_node, access_function);
        //to_console::print_vector(current_dist_vals, 10);
        //std::cout << "\t and directions \t";
        //to_console::print_set(remaining_dirs);
        
        for(auto direction : remaining_dirs)
        {
            
            //std::cout << "Early boundary update: for node " << current_border_node << ", performing entry change " << direction << " -> " << invert_direction(direction);
            //std::cout << " (" << distribution_values[access_function(current_border_node, direction)] << " -> " << distribution_values[access_function(current_border_node, invert_direction(direction))] << " )"<<  std::endl;
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
 * @param velocities a vector containing all velocities
 * @param densities a vector containing all densities
 * @param access_function the access function used to access the distribution values
 */
void bounce_back::perform_inout_boundary_update
(
    border_swap_information &bsi,
    std::vector<double> &distribution_values, 
    // std::vector<velocity> &velocities,
    // std::vector<double> densities,
    access_function access_function
)
{
    std::vector<double> current_dist_vals(DIRECTION_COUNT, 0);
    std::set<unsigned int> remaining_dirs{streaming_directions.begin(), streaming_directions.end()};
    int current_border_node = 0;
    for(auto current : bsi)
    {
        current_border_node = current[0];
        if(std::get<0>(access::get_node_coordinates(current_border_node)) == 1) // Inlet node
        {
            velocity v = INLET_VELOCITY;
            current_dist_vals = maxwell_boltzmann_distribution(v, INLET_DENSITY); // current_dist_vals = maxwell_boltzmann_distribution(velocities[current_border_node], densities[current_border_node]);
            access::set_all_distribution_values
            (
                current_dist_vals,
                distribution_values,
                current_border_node,
                access_function
            );
        }
        else if(std::get<0>(access::get_node_coordinates(current_border_node)) == HORIZONTAL_NODES - 1) // Outlet node
        {
            velocity v = OUTLET_VELOCITY;
            current_dist_vals = maxwell_boltzmann_distribution(v, OUTLET_DENSITY); // current_dist_vals = maxwell_boltzmann_distribution(velocities[current_border_node], densities[current_border_node]);
            access::set_all_distribution_values
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
            //std::cout << "Currently treating node " << current_border_node << std::endl;
            //std::cout << "\t with distribution values \t";
            //current_dist_vals = access::get_distribution_values_of(distribution_values, current_border_node, access_function);
            //to_console::print_vector(current_dist_vals, 10);
            //std::cout << "\t and directions \t";
            //to_console::print_set(remaining_dirs);
            
            for(auto direction : remaining_dirs)
            {
                
                //std::cout << "Early boundary update: for node " << current_border_node << ", performing entry change " << direction << " -> " << invert_direction(direction);
                //std::cout << " (" << distribution_values[access_function(current_border_node, direction)] << " -> " << distribution_values[access_function(current_border_node, invert_direction(direction))] << " )"<<  std::endl;
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
 * @param access_function the access function used to access the distribution values
 */
void bounce_back::update_velocity_input_velocity_output
(
    std::vector<double> &distribution_values, 
    access_function access_function
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
        current_border_node = access::get_node_index(0,y);
        density = macroscopic::density(access::get_distribution_values_of(distribution_values, access::get_neighbor(current_border_node, 5), access_function));
        density = INLET_DENSITY + (INLET_DENSITY - density);
        current_dist_vals = maxwell_boltzmann_distribution(inlet[y-1], density);
        access::set_all_distribution_values
        (
            current_dist_vals,
            distribution_values,
            current_border_node,
            access_function
        );

        // Update outlets
        current_border_node = access::get_node_index(HORIZONTAL_NODES - 1,y);
        density = macroscopic::density(access::get_distribution_values_of(distribution_values, access::get_neighbor(current_border_node, 3), access_function));
        density = OUTLET_DENSITY + (OUTLET_DENSITY - density);
        current_dist_vals = maxwell_boltzmann_distribution(outlet[y-1], density);
        access::set_all_distribution_values
        (
            current_dist_vals,
            distribution_values,
            current_border_node,
            access_function
        );
    }
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
void bounce_back::update_velocity_input_density_output
(
    std::vector<double> &distribution_values, 
    access_function access_function
)
{
    std::vector<double> current_dist_vals(DIRECTION_COUNT, 0);
    int current_border_node = 0;
    velocity v = INLET_VELOCITY;
    double density = 0;

    for(auto y = 0; y < VERTICAL_NODES; ++y)
    {
        // Update inlets
        current_border_node = access::get_node_index(0,y);
        v = INLET_VELOCITY;
        density = INLET_DENSITY;
        current_dist_vals = maxwell_boltzmann_distribution(v, density);
        access::set_all_distribution_values
        (
            current_dist_vals,
            distribution_values,
            current_border_node,
            access_function
        );

        // Update outlets
        current_border_node = access::get_node_index(HORIZONTAL_NODES - 1,y);
        v = macroscopic::flow_velocity(access::get_distribution_values_of(distribution_values, access::get_neighbor(current_border_node, 3), access_function));
        density = OUTLET_DENSITY;
        current_dist_vals = maxwell_boltzmann_distribution(v, density);
        access::set_all_distribution_values
        (
            current_dist_vals,
            distribution_values,
            current_border_node,
            access_function
        );
    }
}

/**
 * @brief Updates the ghost nodes that represent inlet and outlet edges.
 *        When updating, a density border condition will be considered for both the input and the output.
 *        The corresponding values are constants defined in "../include/"defines.hpp".
 * 
 * @param distribution_values a vector containing the distribution values of all nodes
 * @param access_function the access function used to access the distribution values
 */
void bounce_back::update_density_input_density_output
(
    std::vector<double> &distribution_values, 
    access_function access_function
)
{
    std::cout << "Entering" << std::endl;
    std::vector<double> current_dist_vals(DIRECTION_COUNT, 0);
    int current_border_node = 0;
    velocity v = INLET_VELOCITY;
    double density = 0;
    std::cout << "EEverything set up" << std::endl;

    for(auto y = 0; y < VERTICAL_NODES; ++y)
    {
        // Update inlets
        current_border_node = access::get_node_index(0,y);
        v = {0,0};
        density = INLET_DENSITY;
        current_dist_vals = maxwell_boltzmann_distribution(v, density);
        access::set_all_distribution_values
        (
            current_dist_vals,
            distribution_values,
            current_border_node,
            access_function
        );

        // Update outlets
        current_border_node = access::get_node_index(HORIZONTAL_NODES - 1,y);
        v = macroscopic::flow_velocity(access::get_distribution_values_of(distribution_values, access::get_neighbor(current_border_node, 3), access_function));
        density = OUTLET_DENSITY;
        current_dist_vals = maxwell_boltzmann_distribution(v, density);
        access::set_all_distribution_values
        (
            current_dist_vals,
            distribution_values,
            current_border_node,
            access_function
        );
    }
}

/**
 * @brief Initializes all inlet and outlet nodes with their corresponding initial values.
 *        The corresponding values are constants defined in "../include/"defines.hpp".
 * 
 * @param distribution_values a vector containing the distribution values of all nodes
 * @param access_function the access function used to access the distribution values
 */
void bounce_back::initialize_inout
(
    std::vector<double> &distribution_values, 
    access_function access_function
)
{
    std::vector<double> current_dist_vals(DIRECTION_COUNT, 0);
    int current_border_node = 0;
    velocity v = INLET_VELOCITY;
    double density = 0;

    for(auto y = 0; y < VERTICAL_NODES; ++y)
    {
        // Update inlets
        current_border_node = access::get_node_index(0,y);
        v = INLET_VELOCITY;
        density = INLET_DENSITY;
        current_dist_vals = maxwell_boltzmann_distribution(v, density);
        access::set_all_distribution_values
        (
            current_dist_vals,
            distribution_values,
            current_border_node,
            access_function
        );

        // Update outlets
        current_border_node = access::get_node_index(HORIZONTAL_NODES - 1,y);
        v = OUTLET_VELOCITY;
        density = OUTLET_DENSITY;
        current_dist_vals = maxwell_boltzmann_distribution(v, density);
        access::set_all_distribution_values
        (
            current_dist_vals,
            distribution_values,
            current_border_node,
            access_function
        );
    }
}

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