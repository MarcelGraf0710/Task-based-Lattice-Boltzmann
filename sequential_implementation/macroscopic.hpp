#include "defines.hpp"

/**
 * @brief This namespace includes functions for calculating macroscopic observables of non-boundary nodes for the D2Q9I model.
 *        Notice that those functions assume an imcompressible fluid and are false otherwise!
 */
namespace macroscopic
{
    double density(vec_of_dist_val distribution_functions)
    {
        return std::accumulate(distribution_functions.begin(), distribution_functions.end(), 0);
    }

    velocity flow_velocity(vec_of_dist_val distribution_functions)
    {
        double v_x, v_y;
        velocity velocity_vector;

        for(int i = 0; i < DIRECTION_COUNT; ++i)
        {
            velocity_vector = velocity_vectors[i];
            v_x += distribution_functions[i] * velocity_vector[0];
            v_y += distribution_functions[i] * velocity_vector[1];
        }

        velocity v{v_x, v_y};
        return v;
    }
}