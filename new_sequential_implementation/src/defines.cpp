#include "../include/defines.hpp"
#include "../include/utils.hpp"

/** Mapping of directions as proposed by Mattila to the corresponding velocity vectors */
std::map<unsigned int, velocity> velocity_vectors =
{
    {6, {-1, 1}},  {7, {0, 1}}, {8, {1, 1}},   
    {3, {-1, 0}},  {4, {0, 0}},  {5, {1, 0}},   
    {0, {-1, -1}}, {1, {0, -1}}, {2, {1, -1}}   
};

/** Mapping of directions as proposed by Mattila to the weights of the corresponding distribution function */
std::map<unsigned int,double> weights = 
{
    {6, 1.0/36}, {7, 1.0/9}, {8, 1.0/36},
    {3, 1.0/9},  {4, 4.0/9}, {5, 1.0/9},
    {0, 1.0/36}, {1, 1.0/9}, {2, 1.0/36}
};

/** Vector containing the distribution values that actually change within a streaming step */
std::vector<unsigned int> streaming_directions = {0, 1, 2, 3, 5, 6, 7, 8};

/**
 * @brief Returns the Maxwell-Boltzmann-Distribution for all directions in the order proposed by Mattila et al.
 * 
 * @param u two-dimensional velocity vector
 * @param rho density
 * @return the probability of there being a particle with velocity v_direction 
 */
std::vector<double> maxwell_boltzmann_distribution
(
    velocity &u, 
    double rho
)
{
    std::vector<double> result(DIRECTION_COUNT,0);

    for(auto direction = 0; direction < DIRECTION_COUNT; ++direction)
    {
        result[direction] = weights[direction] * rho * 
            (
                1 + 3 * math_utils::dot(velocity_vectors[direction], u) 
                + 9.0/2 * pow(math_utils::dot(velocity_vectors[direction], u), 2)
                - 3.0/2 * math_utils::dot(u, u)
            );
    }
    return result;
}
