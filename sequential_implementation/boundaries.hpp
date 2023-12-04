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
     * 
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
     * @brief Performs the streaming step for a node that borders an inlet within the simulation domain.
     *        Note that in this case, it borders only an inlet and not a wall!
     * 
     * @param density The density at the inlet at this node
     * @param velocity_x The stream velocity at this node in x-direction 
     *                   (note that the velocity in y-direction is assumed to be 0)
     * @param f an array containing all valid incoming distributions
     *          distributions that are not involved may be initialized arbitrarily as they will not be used anyways
     * @return an array containing the correct distribution values for the stream step
     */
    arr_of_dist_val inlet_boundary_stream(arr_of_dist_val f, double velocity_x = INLET_VELOCITY, double density = INLET_DENSITY)
    {
        arr_of_dist_val result = f;
        double rho_times_u_x = 1.0/6 * density * velocity_x;
        result[5] = f[3] + 2.0/3 * density * velocity_x;
        result[8] = f[0] - 0.5 * (f[7] - f[1]) + rho_times_u_x;
        result[2] = f[6] + 0.5 * (f[7] - f[1]) + rho_times_u_x;

        return result;
    }

    /**
     * @brief Performs the streaming step for a node that borders an outlet within the simulation domain.
     *        Note that in this case, it borders only an outlet and not a wall!
     * 
     * @param density The density at the outlet at this node
     * @param velocity_x The stream velocity at this node in x-direction 
     *                   (note that the velocity in y-direction is assumed to be 0)
     * @param f an array containing all valid incoming distributions
     *          distributions that are not involved may be initialized arbitrarily as they will not be used anyways
     * @return an array containing the correct distribution values for the stream step
     */
    arr_of_dist_val outlet_boundary_stream(arr_of_dist_val f, double velocity_x, double density = OUTLET_DENSITY)
    {
        arr_of_dist_val result = f;
        double rho_times_u_x = 1.0/6 * density * velocity_x;
        result[3] = f[5] - 2.0/3 * density * velocity_x;
        result[0] = f[8] + 0.5 * (f[7] - f[1]) - rho_times_u_x;
        result[6] = f[2] - 0.5 * (f[7] - f[1]) - rho_times_u_x;

        return result;
    }

    /**
     * @brief Performs the stream step for a node that is bordered by a lower wall of a simulation domain.
     * 
     * @param wall_velocity this velocity is assumed to be specified at the bordering wall.
     * @param f an array containing all valid incoming distributions
     *          distributions that are not involved may be initialized arbitrarily as they will not be used anyways
     * @return a tuple containing both the correct distribution values and the density (as the latter is needed in calculations anyways)
     */
    std::tuple<arr_of_dist_val, double> lower_wall_boundary_stream(velocity wall_velocity, arr_of_dist_val f)
    {
        arr_of_dist_val result = f;
        double density = 1 / (1 - wall_velocity[1]) * (f[4] + f[5] + f[3] + 2 * (f[1] + f[0] + f[2]));

        double rho_times_u_x = 0.5 * density * wall_velocity[0];
        double rho_times_u_y = 1.0/6 * density * wall_velocity[1];

        result[7] = f[1] + 2.0/3 * density * wall_velocity[1];
        result[8] = f[0] + 0.5 * (f[5] - f[3]) + rho_times_u_x + rho_times_u_y;
        result[6] = f[2] + 0.5 * (f[5] - f[3]) - rho_times_u_x + rho_times_u_y;

        std::tuple<arr_of_dist_val, double> result_tuple(result, density);
        return result_tuple;
    }

    /**
     * @brief Performs the stream step for a node that is bordered by an upper wall of a simulation domain.
     * 
     * @param wall_velocity this velocity is assumed to be specified at the bordering wall.
     * @param f an array containing all valid incoming distributions
     *          distributions that are not involved may be initialized arbitrarily as they will not be used anyways
     * @return a tuple containing both the correct distribution values and the density (as the latter is needed in calculations anyways)
     */
    std::tuple<arr_of_dist_val, double> upper_wall_boundary_stream(velocity wall_velocity, arr_of_dist_val f)
    {
        arr_of_dist_val result = f;
        double density = 1 / (1 - wall_velocity[1]) * (f[4] + f[5] + f[3] + 2 * (f[1] + f[0] + f[2]));

        double rho_times_u_x = 0.5 * density * wall_velocity[0];
        double rho_times_u_y = 1.0/6 * density * wall_velocity[1];

        result[1] = f[7] - 2.0/3 * density * wall_velocity[1];
        result[0] = f[8] - 0.5 * (f[3] - f[5]) - rho_times_u_x - rho_times_u_y;
        result[2] = f[6] + 0.5 * (f[3] - f[5]) + rho_times_u_x - rho_times_u_y;

        std::tuple<arr_of_dist_val, double> result_tuple(result, density);
        return result_tuple;
    }

    /**
     * @brief Performs the stream step for a node that is bordered by a lower wall and an inlet of a simulation domain.
     * 
     * @param f an array containing all valid incoming distributions
     *          distributions that are not involved may be initialized arbitrarily as they will not be used anyways
     * @param density the density at the inlet
     * @return an array containing the correct distribution values for the stream step
     */
    arr_of_dist_val lower_inlet_boundary_stream(arr_of_dist_val f, double density = INLET_DENSITY)
    {
        arr_of_dist_val result = f;
        result[5] = f[3];
        result[7] = f[1];
        result[8] = f[0];
        result[2] = 0.5 * (density - (result[4] + 2 * result[0] + 2 * result[1] + 2 * result[3]));
        result[6] = f[2];
        return result;
    }

    /**
     * @brief Performs the stream step for a node that is bordered by an upper wall and an inlet of a simulation domain.
     * 
     * @param f an array containing all valid incoming distributions
     *          distributions that are not involved may be initialized arbitrarily as they will not be used anyways
     * @param density the density at the inlet
     * @return an array containing the correct distribution values for the stream step
     */
    arr_of_dist_val upper_inlet_boundary_stream(arr_of_dist_val f, double density = INLET_DENSITY)
    {
        arr_of_dist_val result = f;
        result[5] = f[3];
        result[1] = f[7];
        result[2] = f[6];
        result[0] = 0.5 * (density - (result[4] + 2 * result[3] + 2 * result[6] + 2 * result[7]));
        result[8] = f[0];
        return result;
    }

    /**
     * @brief Performs the stream step for a node that is bordered by a lower wall and an outlet of a simulation domain.
     * 
     * @param f an array containing all valid incoming distributions
     *          distributions that are not involved may be initialized arbitrarily as they will not be used anyways
     * @param density the density at the outlet
     * @return an array containing the correct distribution values for the stream step
     */
    arr_of_dist_val lower_outlet_boundary_stream(arr_of_dist_val f, double density = INLET_DENSITY)
    {
        arr_of_dist_val result = f;
        result[3] = f[5];
        result[7] = f[1];
        result[6] = f[2];
        result[0] = 0.5 * (density - (result[4] + 2 * result[1] + 2 * result[2] + 2 * result[5]));
        result[8] = f[0];
        return result;
    }

    /**
     * @brief Performs the stream step for a node that is bordered by an upper wall and an outlet of a simulation domain.
     * 
     * @param f an array containing all valid incoming distributions
     *          distributions that are not involved may be initialized arbitrarily as they will not be used anyways
     * @param density the density at the outlet
     * @return an array containing the correct distribution values for the stream step
     */
    arr_of_dist_val upper_outlet_boundary_stream(arr_of_dist_val f, double density = INLET_DENSITY)
    {
        arr_of_dist_val result = f;
        result[3] = f[5];
        result[1] = f[7];
        result[0] = f[8];
        result[2] = 0.5 * (density - (result[4] + 2 * result[5] + 2 * result[7] + 2 * result[8]));
        result[6] = f[2];
        return result;
    }

    /**
     * @brief Returns the according result after a streaming step for the specified boundary scenario.
     *        So far, this method is only conceptional and it is likely to be removed in the future for the sake of complexity reduction.
     * 
     * @param boundary_spec any tuple included within the boundary_scenarios namespace               
     * @return a tuple containing an array with all distribution function values and the density value
     */
    std::tuple<arr_of_dist_val, double> perform_boundary_stream(boundary_tuple boundary_spec, velocity wall_velocity, arr_of_dist_val f)
    {
        if(boundary_spec == boundary_scenarios::inlet)
        {
            throw std::runtime_error(std::string("This boundary scenario is not implemented yet!"));
        }
        else if(boundary_spec == boundary_scenarios::outlet)
        {
            throw std::runtime_error(std::string("This boundary scenario is not implemented yet!"));
        }
        else if(boundary_spec == boundary_scenarios::wall_up)
        {
            throw std::runtime_error(std::string("This boundary scenario is not implemented yet!"));
        }
        else if(boundary_spec == boundary_scenarios::wall_down)
        {
            throw std::runtime_error(std::string("This boundary scenario is not implemented yet!"));
        }
        else if(boundary_spec == boundary_scenarios::lower_inlet)
        {
            throw std::runtime_error(std::string("This boundary scenario is not implemented yet!"));
        }
        else if(boundary_spec == boundary_scenarios::upper_inlet)
        {
            throw std::runtime_error(std::string("This boundary scenario is not implemented yet!"));
        }
        else if(boundary_spec == boundary_scenarios::lower_outlet)
        {
            throw std::runtime_error(std::string("This boundary scenario is not implemented yet!"));
        }
        else if(boundary_spec == boundary_scenarios::upper_outlet)
        {
           throw std::runtime_error(std::string("This boundary scenario is not implemented yet!")); 
        }
        else
        {
            throw std::runtime_error(std::string("This boundary scenario is not supported!"));
        }
    }
}