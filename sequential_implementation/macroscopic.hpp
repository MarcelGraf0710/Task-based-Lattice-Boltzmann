#include "defines.hpp"

/**
 * @brief This namespace includes functions for calculating macroscopic observables of non-boundary nodes for the D2Q9I model.
 *        Notice that those functions assume an imcompressible fluid and are false otherwise!
 */
namespace macroscopic
{
    double density(arr_of_dist_val distribution_functions)
    {
        return std::accumulate(distribution_functions.begin(), distribution_functions.end(), 0);
    }

    velocity flow_velocity(arr_of_dist_val distribution_functions)
    {
        double v_x, v_y;

        for(auto velocity_vector : velocity_vectors)
        {
            v_x += velocity_vector.second[0];
            v_y += velocity_vector.second[1];
        }

        velocity v{v_x, v_y};
        return v;
    }
}