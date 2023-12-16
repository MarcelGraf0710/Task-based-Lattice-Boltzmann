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

std::vector<unsigned int> streaming_directions = {0, 1, 2, 3, 5, 6, 7, 8};

std::vector<double> maxwell_boltzmann_distribution(const velocity &u, double rho)
{
    std::vector<double> result;
    for(auto direction = 0; direction < DIRECTION_COUNT; ++direction)
    {
        result.push_back(weights[direction] * rho * (1 + 3 * math_utils::dot(velocity_vectors[direction], u) 
                    + 9.0/2 * pow(math_utils::dot(velocity_vectors[direction], u), 2)
                    - 3.0/2 * math_utils::dot(velocity_vectors[direction], u)));
    }
    return result;
}
