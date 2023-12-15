#include "../include/defines.hpp"
#include "../include/utils.hpp"
#include <iostream>

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

double maxwell_boltzmann_distribution(velocity &u, double rho, unsigned int direction)
{
    return weights[direction] * rho * (1 + 3 * math_utils::dot(velocity_vectors[direction], u) 
    + 9.0/2 * pow(math_utils::dot(velocity_vectors[direction], u), 2) - 3.0/2 * math_utils::dot(velocity_vectors[direction], u));
}

vec_of_dist_val maxwell_boltzmann_distribution(velocity &u, double rho)
{
    std::cout << "Entering maxwell_boltzmann_distribution" << std::endl;
    vec_of_dist_val result;
    for(auto direction = 0; direction < DIRECTION_COUNT; ++direction)
    {
        std::cout << "weights[" << direction << "] = " << weights[direction] << std::endl;
        std::cout << "velocity_vectors[" << direction << "] = " << "(" << velocity_vectors[direction][0] << ", " << velocity_vectors[direction][1] << ")" << std::endl;
        result.push_back(weights[direction] * rho * (1 + 3 * math_utils::dot(velocity_vectors[direction], u) 
                    + 9.0/2 * pow(math_utils::dot(velocity_vectors[direction], u), 2)
                    - 3.0/2 * math_utils::dot(velocity_vectors[direction], u)));
    }
    std::cout << "Leaving maxwell_boltzmann_distribution" << std::endl;
    return result;
}
