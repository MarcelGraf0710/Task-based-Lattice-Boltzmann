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

std::vector<double> maxwell_boltzmann_distribution(velocity &u, double rho)
{
    //std::cout << "Entering mb distro with params u = (" << u[0] << ", " << u[1] << "), rho = " << rho << std::endl;
    std::vector<double> result;
    result.reserve(DIRECTION_COUNT);
    for(auto direction = 0; direction < DIRECTION_COUNT; ++direction)
    {
        //std::cout << "mb distro @ dir " << direction << std::endl;
        //std::cout << "weights[direction] = " << weights[direction] << std::endl;
        //std::cout << "velocity_vectors[direction] = (" << velocity_vectors[direction][0] << ", " << velocity_vectors[direction][1] << ")" << std::endl;
        //std::cout << "Dot product = " << math_utils::dot(velocity_vectors[direction], u) << std::endl;
        //std::cout << "Length of result is " << result.size() << std::endl;
        //result.push_back(666);
        //std::cout << "Got: " << *result.end() << std::endl;
        result.push_back(weights[direction] * rho * (1 + 3 * math_utils::dot(velocity_vectors[direction], u) 
                    + 9.0/2 * pow(math_utils::dot(velocity_vectors[direction], u), 2)
                    - 3.0/2 * math_utils::dot(velocity_vectors[direction], u)));
    }
    return result;
}
