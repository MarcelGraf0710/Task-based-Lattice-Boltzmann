/* Include directives */
#include <vector>
#include <list>
#include <complex>
#include <map>
#include "utils.hpp"
#include <functional>
#include <valarray>
#include <numeric>

/* General definitions (NOTE: this may be outsourced to a .csv file in the future to increase modifiability.) */
#define DIMENSION_COUNT 2
#define DIRECTION_COUNT 9
#define TIME_STEP 1
#define SPACE_STEP 1
#define VERTICAL_NODES 3
#define HORIZONTAL_NODES 5
#define TOTAL_NODE_COUNT VERTICAL_NODES * HORIZONTAL_NODES
#define SPEED_OF_SOUND 1.0/9
#define TEMPERATURE 20
#define BOLTZMANN_CONSTANT 1.380649e-23
#define RELAXATION_TIME 0.8

/* Convenience and readability type definition */
typedef std::array<std::array<double, DIMENSION_COUNT>, DIRECTION_COUNT> arr_of_v; // Array of velocity vectors
typedef std::array<double, DIMENSION_COUNT> velocity; // velocity vector
typedef std::vector<double> vec_of_dist_val; // array of distribution function values
typedef std::vector<double> all_distributions; // A vector containing all distribution values

/**
 * @brief This type represents a tuple containing all boundary nodes in the following order:
 *        Corner nodes (upper left, lower left, upper right, lower right) ->
 *        inlet nodes ->
 *        outlet nodes ->
 *        upper wall nodes ->
 *        lower wall nodes.
 *        In order to avoid unnecessary branching for comparably few nodes, these nodes will be treated separately by design.
 */
typedef std::tuple<
                    std::array<int, 4>, 
                    std::array<int, VERTICAL_NODES - 2>,
                    std::array<int, VERTICAL_NODES - 2>,
                    std::array<int, HORIZONTAL_NODES - 2>,
                    std::array<int, HORIZONTAL_NODES - 2>
                    > 
                node_tuple;

/**
 * @brief This type stands for an access function. Node values can be stored in different layout and via this function,
 *        the corresponding access scheme can be specified.
 */
typedef std::function<unsigned int(unsigned int, unsigned int)> access_function;

/* Inlet and outlet behavior */
#define INLET_VELOCITY 1.0/100
#define INLET_DENSITY 2
#define OUTLET_DENSITY 2 

/* Velocity vectors */
std::map<int, velocity> velocity_vectors 
{
    {6, {-1, 1}},  {7, {0, -1}}, {8, {1, 1}},   
    {3, {-1, 0}},  {4, {0, 0}},  {5, {1, 0}},   
    {0, {-1, -1}}, {1, {0, -1}}, {2, {1, -1}}   
};

/* Weights */
std::map<int,double> weights
    {
        {6, 1.0/36}, {7, 1.0/9}, {8, 1.0/36},
        {3, 1.0/9},  {4, 4.0/9}, {5, 1.0/9},
        {0, 1.0/36}, {1, 1.0/9}, {2, 1.0/36}
    };

/**
 * @brief The Maxwell-Boltzmann-Distribution marks the equilibrium distribution of particles.
 * 
 * @param u two-dimensional velocity vector
 * @param rho density
 * @param direction direction according to the scheme proposed by Mattila et al.
 * @return the probability of there being a particle with velocity v_direction 
 */
double maxwell_boltzmann_distribution(velocity u, double rho, unsigned int direction)
{
    return weights[direction] * rho * (1 + 3 * math_utils::dot(velocity_vectors[direction], u) 
    + 9.0/2 * pow(math_utils::dot(velocity_vectors[direction], u), 2) - 3.0/2 * math_utils::dot(velocity_vectors[direction], u));
}

/**
 * @brief Returns the Maxwell-Boltzmann-Distribution for all directions in the order proposed by Mattila et al.
 * 
 * @param u two-dimensional velocity vector
 * @param rho density
 * @return the probability of there being a particle with velocity v_direction 
 */
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
