#ifndef BOUNDARIES_HPP
#define BOUNDARIES_HPP

#include "access.hpp"
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
        {boundary_scenarios::inlet, {1,2,5,7,8}},
        {boundary_scenarios::outlet, {0,1,3,6,7}},
        {boundary_scenarios::lower_inlet, {5,7,8}},
        {boundary_scenarios::upper_inlet, {1,2,5}},
        {boundary_scenarios::lower_outlet, {3,6,7}},
        {boundary_scenarios::upper_outlet, {0,1,3}},
        {boundary_scenarios::wall_up, {0,1,2,3,5}},
        {boundary_scenarios::wall_down, {3,5,6,7,8}},
    };

    /**
     * @brief Performs the after-streaming value update for a node that borders an inlet within the simulation domain.
     *        Note that in this case, it borders only an inlet and not a wall!
     *      
     * @param destination the vector containing the distribution values for the current node
     * @param access_function the access pattern of the corresponding destination vector
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
        double density = INLET_DENSITY);

    /**
     * @brief Performs the after-streaming value update for a node that borders an outlet within the simulation domain.
     *        Note that in this case, it borders only an outlet and not a wall!
     *      
     * @param destination the vector containing the distribution values for the current node
     * @param access_function the access pattern of the corresponding destination vector
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
        double density = OUTLET_DENSITY);

    /**
     * @brief Performs the after-streaming value update for a node that borders a lower wall within the simulation domain.
     *        Note that in this case, it borders only a lower wall and not an inlet or outlet!
     *      
     * @param destination the vector containing the distribution values for the current node
     * @param access_function the access pattern of the corresponding destination vector
     * @param x the x coordinate of the node (given that the y coordinate is 0)
     * @param wall_velocity the velocity at the wall (default is {0,0})
     * @return the density at this node (as it is needed for calculations anyways)
     */
    double lower_wall_boundary_stream(  
        all_distributions &source,      
        all_distributions &destination, 
        access_function &access_function, 
        int x,
        velocity wall_velocity = {0,0});

    /**
     * @brief Performs the after-streaming value update for a node that borders an upper wall within the simulation domain.
     *        Note that in this case, it borders only an upper wall and not an inlet or outlet!
     *      
     * @param destination the vector containing the distribution values for the current node
     * @param access_function the access pattern of the corresponding destination vector
     * @param x the x coordinate of the node (given that the y coordinate is VERTICAL_NODES - 1)
     * @param wall_velocity the velocity at the wall (default is {0,0})
     * @return the density at this node (as it is needed for calculations anyways)
     */
    double upper_wall_boundary_stream(   
        all_distributions &source,      
        all_distributions &destination, 
        access_function &access_function, 
        int x, 
        velocity wall_velocity = {0,0});

    /**
     * @brief Performs the after-stream update for a node that is bordered by a lower wall and an inlet of a simulation domain.
     * 
     * @param access_function the access pattern of the corresponding destination vector
     * @param destination the vector containing the distribution values for the current node
     * @param density the density at this node (default value is INLET_DENSITY)
     */
    void lower_inlet_boundary_stream(
        all_distributions &source, 
        all_distributions &destination,
        access_function access_function,
        double density = INLET_DENSITY);

    /**
     * @brief Performs the after-stream update for a node that is bordered by an upper wall and an inlet of a simulation domain.
     * 
     * @param access_function the access pattern of the corresponding destination vector
     * @param destination the vector containing the distribution values for the current node
     * @param density the density at this node (default value is INLET_DENSITY)
     */
    void upper_inlet_boundary_stream(
        all_distributions &source, 
        all_distributions &destination, 
        access_function access_function,
        double density = INLET_DENSITY);

    /**
     * @brief Performs the after-stream update for a node that is bordered by a lower wall and an outlet of a simulation domain.
     * 
     * @param access_function the access pattern of the corresponding destination vector
     * @param destination the vector containing the distribution values for the current node
     * @param density the density at this node (default value is INLET_DENSITY as the general assumption is an uncompressible stream)
     */
    void lower_outlet_boundary_stream(
        all_distributions &source,        
        all_distributions &destination,        
        access_function access_function,
        double density = INLET_DENSITY);

    /**
     * @brief Performs the after-stream update for a node that is bordered by an upper wall and an outlet of a simulation domain.
     * 
     * @param access_function the access pattern of the corresponding destination vector
     * @param destination the vector containing the distribution values for the current node
     * @param density the density at this node (default value is INLET_DENSITY as the general assumption is an uncompressible stream)
     */
    void upper_outlet_boundary_stream(   
        all_distributions &source, 
        all_distributions &destination,         
        access_function access_function,
        double density = INLET_DENSITY);
}

#endif