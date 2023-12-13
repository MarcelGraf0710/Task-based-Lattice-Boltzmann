#include "../include/defines.hpp"
#include "../include/utils.hpp"

/* Velocity vectors */
std::map<int, velocity> velocity_vectors =
{
    {6, {-1, 1}},  {7, {0, -1}}, {8, {1, 1}},   
    {3, {-1, 0}},  {4, {0, 0}},  {5, {1, 0}},   
    {0, {-1, -1}}, {1, {0, -1}}, {2, {1, -1}}   
};

/* Weights */
std::map<int,double> weights = 
    {
        {6, 1.0/36}, {7, 1.0/9}, {8, 1.0/36},
        {3, 1.0/9},  {4, 4.0/9}, {5, 1.0/9},
        {0, 1.0/36}, {1, 1.0/9}, {2, 1.0/36}
    };

double maxwell_boltzmann_distribution(velocity u, double rho, unsigned int direction)
{
    return weights[direction] * rho * (1 + 3 * math_utils::dot(velocity_vectors[direction], u) 
    + 9.0/2 * pow(math_utils::dot(velocity_vectors[direction], u), 2) - 3.0/2 * math_utils::dot(velocity_vectors[direction], u));
}

vec_of_dist_val maxwell_boltzmann_distribution(velocity u, double rho)
{
    vec_of_dist_val result;
    for(auto direction = 0; direction < DIRECTION_COUNT; ++direction)
    {
        result[direction] = weights[direction] * rho * (1 + 3 * math_utils::dot(velocity_vectors[direction], u) 
                    + 9.0/2 * pow(math_utils::dot(velocity_vectors[direction], u), 2)
                    - 3.0/2 * math_utils::dot(velocity_vectors[direction], u));
    }
    return result;
}
