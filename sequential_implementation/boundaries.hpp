#include "defines.hpp"
#include <optional>

/**
 * @brief This namespace contains all function representations of boundary conditions used in the lattice-Boltzmann model.
 */
namespace boundaries
{
    /**
     * @brief Convenience enum that allows for easier and more readable boundary condition differenciation.
     *        The intended use is with an initializer list that filters out the respective boundary situation.
     */
    enum boundary_type {INLET, OUTLET, WALL, VELOCITY, PRESSURE};
    enum boundary_direction {UP, DOWN, LEFT, RIGHT};
    typedef std::tuple<std::tuple<boundary_type, boundary_direction>, std::optional<std::tuple<boundary_type, boundary_direction>>> boundary_tuple;

    /**
     * @brief This namespace contains all currently implemented boundary node scenarios. 
     *        They are used in the determination process of the correct boundary treatment.
     */
    namespace boundary_scenarios
    {
        const boundary_tuple inlet {std::tuple(INLET, LEFT), std::nullopt};
        const boundary_tuple outlet {std::tuple(OUTLET, RIGHT), std::nullopt};
        const boundary_tuple lower_inlet {std::tuple(INLET, LEFT), std::tuple(WALL, DOWN)};
        const boundary_tuple upper_inlet {std::tuple(INLET, LEFT), std::tuple(WALL, UP)};
        const boundary_tuple lower_outlet {std::tuple(INLET, LEFT), std::tuple(WALL, DOWN)};
        const boundary_tuple upper_outlet {std::tuple(OUTLET, RIGHT), std::tuple(WALL, UP)};
        const boundary_tuple wall_up {std::tuple(WALL, UP), std::nullopt};
        const boundary_tuple wall_down {std::tuple(WALL, DOWN), std::nullopt};
    };

    /**
     * @brief This map contains the directions that each kind of border node is facing a neighbor node,
     *        i.e. it contains all directions that perform actual streams.
     */
    std::map<boundary_tuple, std::list<int>> neighbor_directions
    {
        {inlet, {1,2,5,7,8}},
        {outlet, {0,1,3,6,7}},
        {lower_inlet, {5,7,8}},
        {upper_inlet, {1,2,5}},
        {lower_outlet, {3,6,7}},
        {upper_outlet, {0,1,3}},
        {wall_up, {0,1,2,3,5}},
        {wall_down, {3,5,6,7,8}},
    };

    /**
     * @brief Utility function that gets all distribution values from the specified source using the specified access function.
     * 
     * @param source this array contains the distribution values of all nodes
     * @param access_function this access function is used to retrieve the distribution values from the source array
     * @param node_index this is the index the node has within the access pattern
     * @return an array containing all distribution values
     */
    inline arr_of_dist_val get_all_distribution_values(
        all_distributions &source, 
        access_function &access_function, 
        unsigned int node_index)
    {
        arr_of_dist_val result;
        for (int direction = 0; direction < DIRECTION_COUNT; ++direction)
        {
            result[direction] = source[access_function(node_index, direction)];
        }
        return result;
    }

    /**
     * @brief Utility function to update all specified distribution values at the specified destination using the specified access function.
     * 
     * @tparam d defines the length of the array of changes (should be assigned automatically, can be ignored)
     * @param changes an array containing all directions that were changed
     * @param destination this array contains the distribution values of all nodes
     * @param access_function this access function is used to determine the indices of the directions within the destination array
     * @param node_index the index of the node whose distribution values are updated
     * @param dist_values an array containing all distribution values of the specified node (can be optimized, I know...)
     */
    template <unsigned long d>
    inline void write_updated_border_values(
        std::array<unsigned int, d> &changes, 
        all_distributions &destination, 
        access_function &access_function, 
        unsigned int node_index, 
        arr_of_dist_val &dist_values)
    {
        for (auto direction : changes)
        {
            destination[access_function(node_index, direction)] = dist_values[direction];
        }
    }

    /**
     * @brief Performs the after-streaming value update for a node that borders an inlet within the simulation domain.
     *        Note that in this case, it borders only an inlet and not a wall!
     *      
     * @param destination the array containing the distribution values for the current node
     * @param access_function the access pattern of the corresponding destination array
     * @param y the y coordinate of the inlet node (given that the x coordinate is 0)
     * @param velocity_x the velocity in x direction of this node (default value is INLET_VELOCITY)
     * @param density the density at this node (default value is INLET_DENSITY)
     */
    void inlet_boundary_stream(
        all_distributions &source, 
        all_distributions &destination, 
        access_function &access_function, 
        int y,
        double velocity_x = INLET_VELOCITY, 
        double density = INLET_DENSITY)
    {
        arr_of_dist_val result;
        unsigned int node_index = access::get_node_index(0, y);
        result = get_all_distribution_values(source, access_function, node_index);
        double rho_times_u_x = 1.0/6 * density * velocity_x;

        std::array<unsigned int, 3> changes = {5, 8, 2};
        result[5] = result[3] + 2.0/3 * density * velocity_x;
        result[8] = result[0] - 0.5 * (result[7] - result[1]) + rho_times_u_x;
        result[2] = result[6] + 0.5 * (result[7] - result[1]) + rho_times_u_x;
        
        write_updated_border_values(changes, destination, access_function, node_index, result);
    }

    /**
     * @brief Performs the after-streaming value update for a node that borders an outlet within the simulation domain.
     *        Note that in this case, it borders only an outlet and not a wall!
     *      
     * @param destination the array containing the distribution values for the current node
     * @param access_function the access pattern of the corresponding destination array
     * @param y the y coordinate of the outlet node (given that the x coordinate is HORIZONTAL_NODES - 1)
     * @param velocity_x the velocity in x direction of this node
     * @param density the density at this node (default value is INLET_DENSITY as the general assumption is an uncompressible stream)
     */
    void outlet_boundary_stream(
        all_distributions &source, 
        all_distributions &destination, 
        access_function &access_function, 
        int y, 
        double velocity_x, 
        double density = OUTLET_DENSITY)
    {
        arr_of_dist_val result;
        unsigned int node_index = access::get_node_index(HORIZONTAL_NODES - 1, y);
        result = get_all_distribution_values(source, access_function, node_index);
        double rho_times_u_x = 1.0/6 * density * velocity_x;

        std::array<unsigned int, 3> changes = {5, 8, 2};
        result[3] = result[5] - 2.0/3 * density * velocity_x;
        result[0] = result[8] + 0.5 * (result[7] - result[1]) - rho_times_u_x;
        result[6] = result[2] - 0.5 * (result[7] - result[1]) - rho_times_u_x;

        write_updated_border_values(changes, destination, access_function, node_index, result);
    }

    /**
     * @brief Performs the after-streaming value update for a node that borders a lower wall within the simulation domain.
     *        Note that in this case, it borders only a lower wall and not an inlet or outlet!
     *      
     * @param destination the array containing the distribution values for the current node
     * @param access_function the access pattern of the corresponding destination array
     * @param x the x coordinate of the node (given that the y coordinate is 0)
     * @param wall_velocity the velocity at the wall (default is {0,0})
     * @return the density at this node (as it is needed for calculations anyways)
     */
    double lower_wall_boundary_stream(  
        all_distributions &source,      
        all_distributions &destination, 
        access_function &access_function, 
        int x, 
        velocity wall_velocity = {0,0})
    {
        arr_of_dist_val result;
        unsigned int node_index = access::get_node_index(x, 0);
        result = get_all_distribution_values(source, access_function, node_index);
        double density = 1 / (1 - wall_velocity[1]) * (result[4] + result[5] + result[3] + 2 * (result[1] + result[0] + result[2]));

        double rho_times_u_x = 0.5 * density * wall_velocity[0];
        double rho_times_u_y = 1.0/6 * density * wall_velocity[1];

        std::array<unsigned int, 3> changes = {6,7,8};
        result[7] = result[1] + 2.0/3 * density * wall_velocity[1];
        result[8] = result[0] + 0.5 * (result[5] - result[3]) + rho_times_u_x + rho_times_u_y;
        result[6] = result[2] + 0.5 * (result[5] - result[3]) - rho_times_u_x + rho_times_u_y;

        write_updated_border_values(changes, destination, access_function, node_index, result);
        return density;
    }

    /**
     * @brief Performs the after-streaming value update for a node that borders an upper wall within the simulation domain.
     *        Note that in this case, it borders only an upper wall and not an inlet or outlet!
     *      
     * @param destination the array containing the distribution values for the current node
     * @param access_function the access pattern of the corresponding destination array
     * @param x the x coordinate of the node (given that the y coordinate is VERTICAL_NODES - 1)
     * @param wall_velocity the velocity at the wall (default is {0,0})
     * @return the density at this node (as it is needed for calculations anyways)
     */
    double upper_wall_boundary_stream(   
        all_distributions &source,      
        all_distributions &destination, 
        access_function &access_function, 
        int x, 
        velocity wall_velocity = {0,0})
    {
        arr_of_dist_val result;
        unsigned int node_index = access::get_node_index(x, VERTICAL_NODES - 1);
        result = get_all_distribution_values(source, access_function, node_index);
        double density = 1 / (1 - wall_velocity[1]) * (result[4] + result[5] + result[3] + 2 * (result[1] + result[0] + result[2]));
        double density = 1 / (1 - wall_velocity[1]) * (result[4] + result[5] + result[3] + 2 * (result[1] + result[0] + result[2]));

        double rho_times_u_x = 0.5 * density * wall_velocity[0];
        double rho_times_u_y = 1.0/6 * density * wall_velocity[1];

        std::array<unsigned int, 3> changes = {0,1,2};
        result[1] = result[7] - 2.0/3 * density * wall_velocity[1];
        result[0] = result[8] - 0.5 * (result[3] - result[5]) - rho_times_u_x - rho_times_u_y;
        result[2] = result[6] + 0.5 * (result[3] - result[5]) + rho_times_u_x - rho_times_u_y;

        write_updated_border_values(changes, destination, access_function, node_index, result);
        return density;
    }

    /**
     * @brief Performs the after-stream update for a node that is bordered by a lower wall and an inlet of a simulation domain.
     * 
     * @param access_function the access pattern of the corresponding destination array
     * @param destination the array containing the distribution values for the current node
     * @param density the density at this node (default value is INLET_DENSITY)
     */
    void lower_inlet_boundary_stream(
        all_distributions &source, 
        all_distributions &destination,
        access_function access_function,
        double density = INLET_DENSITY)
    {
        arr_of_dist_val result;
        unsigned int lower_inlet_node = access::get_node_index(0, 0);
        result = get_all_distribution_values(source, access_function, lower_inlet_node);

        std::array<unsigned int, 5> changes = {5, 7, 8, 2, 6};
        result[5] = result[3];
        result[7] = result[1];
        result[8] = result[0];
        result[2] = 0.5 * (density - (result[4] + 2 * result[0] + 2 * result[1] + 2 * result[3]));
        result[6] = result[2];

        write_updated_border_values(changes, destination, access_function, lower_inlet_node, result);
    }

    /**
     * @brief Performs the after-stream update for a node that is bordered by an upper wall and an inlet of a simulation domain.
     * 
     * @param access_function the access pattern of the corresponding destination array
     * @param destination the array containing the distribution values for the current node
     * @param density the density at this node (default value is INLET_DENSITY)
     */
    void upper_inlet_boundary_stream(
        all_distributions &source, 
        all_distributions &destination, 
        access_function access_function,
        double density = INLET_DENSITY)
    {
        arr_of_dist_val result;
        unsigned int upper_inlet_node = access::get_node_index(0, HORIZONTAL_NODES - 1);
        result = get_all_distribution_values(source, access_function, upper_inlet_node);

        std::array<unsigned int, 5> changes = {5, 1, 2, 0, 8};
        result[5] = result[3];
        result[1] = result[7];
        result[2] = result[6];
        result[0] = 0.5 * (density - (result[4] + 2 * result[3] + 2 * result[6] + 2 * result[7]));
        result[8] = result[0];

        write_updated_border_values(changes, destination, access_function, upper_inlet_node, result);
    }

    /**
     * @brief Performs the after-stream update for a node that is bordered by a lower wall and an outlet of a simulation domain.
     * 
     * @param access_function the access pattern of the corresponding destination array
     * @param destination the array containing the distribution values for the current node
     * @param density the density at this node (default value is INLET_DENSITY as the general assumption is an uncompressible stream)
     */
    void lower_outlet_boundary_stream(
        all_distributions &source,        
        all_distributions &destination,        
        access_function access_function,

        double density = INLET_DENSITY)
    {
        arr_of_dist_val result;
        unsigned int lower_outlet_node = access::get_node_index(VERTICAL_NODES - 1, 0);
        result = get_all_distribution_values(source, access_function, lower_outlet_node);
        
        std::array<unsigned int, 5> changes = {3,7,6,0,8};
        result[3] = result[5];
        result[7] = result[1];
        result[6] = result[2];
        result[0] = 0.5 * (density - (result[4] + 2 * result[1] + 2 * result[2] + 2 * result[5]));
        result[8] = result[0];

        write_updated_border_values(changes, destination, access_function, lower_outlet_node, result);
    }

    /**
     * @brief Performs the after-stream update for a node that is bordered by an upper wall and an outlet of a simulation domain.
     * 
     * @param access_function the access pattern of the corresponding destination array
     * @param destination the array containing the distribution values for the current node
     * @param density the density at this node (default value is INLET_DENSITY as the general assumption is an uncompressible stream)
     */
    void upper_outlet_boundary_stream(   
        all_distributions &source, 
        all_distributions &destination,         
        access_function access_function,
        double density = INLET_DENSITY)
    {
        arr_of_dist_val result;
        unsigned int lower_outlet_node = access::get_node_index(VERTICAL_NODES - 1, HORIZONTAL_NODES - 1);
        result = get_all_distribution_values(source, access_function, lower_outlet_node);
        
        std::array<unsigned int, 5> changes = {3,1,0,2,6};
        result[3] = result[5];
        result[1] = result[7];
        result[0] = result[8];
        result[2] = 0.5 * (density - (result[4] + 2 * result[5] + 2 * result[7] + 2 * result[8]));
        result[6] = result[2];

        write_updated_border_values(changes, destination, access_function, lower_outlet_node, result);
    }
}