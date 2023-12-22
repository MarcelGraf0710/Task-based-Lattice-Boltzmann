#ifndef DEFINES_HPP
#define DEFINES_HPP

/* Include directives */

#include <vector>
#include <list>
#include <complex>
#include <map>
#include <functional>
#include <valarray>
#include <numeric>

/* General definitions (NOTE: this may be outsourced to a .csv file in the future to increase modifiability.) */

#define DIMENSION_COUNT 2
#define DIRECTION_COUNT 9
#define TIME_STEP 1
#define SPACE_STEP 1
#define VERTICAL_NODES 5
#define HORIZONTAL_NODES 7
#define TOTAL_NODE_COUNT VERTICAL_NODES * HORIZONTAL_NODES
#define SPEED_OF_SOUND 1.0/9
#define BOLTZMANN_CONSTANT 1.380649e-23
#define RELAXATION_TIME 0.8

/* Convenience and readability type definitions */

/**
 * @brief Representation of a velocity vector
 */
typedef std::array<double, DIMENSION_COUNT> velocity; 

/**
 * @brief Convenience type definition that represents a vector from which the boundary treatment 
 *        of all nodes can be retrieved. Information is stored in the following way:
 * 
 *        Each entry of the outer vector represents a border node, called BORDER_NODE for explanation.
 *        Every BORDER_NODE is a vector of tuples t_i = (NEIGHBOR_i, DIRECTION_i).
 *        Each t_i represents an adjacency relation of BORDER_NODE to a ghost node (including solid nodes).
 *        The first tuple, t_0, is always a self-pointing entry (BORDER_NODE, 4) from which the index 
 *        of the border node is restored.
 *        Any further tuples t_1, ..., t_x store the information what nodes BORDER_NODE is adjacent to 
 *        and by which velocity vector the respective neighbor can be reached.
 * 
 *        This is used for the halfway bounce-back boundary treatment where BORDER_NODE will copy the 
 *        respective f_(DIRECTION_i) from NEIGHBOR_i.
 */
typedef std::vector<std::vector<std::tuple<unsigned int, unsigned int>>> border_adjacency;

/**
 * @brief Convenience type definition that represents a vector from which the boundary treatment 
 *        of all nodes can be retrieved. Information is stored in the following way:
 *        Each entry of the outer vector represents a border node, called BORDER_NODE for explanation.
 *        Every BORDER_NODE is a vector with the following information:
 *        - 0th entry: The index of BORDER_NODE
 *        - Further entries: The directions pointing to ghost nodes (including solid nodes).
 * 
 *        This is used for the halfway bounce-back boundary treatment where BORDER_NODE will copy the 
 *        respective distribution values BEFORE streaming.
 */
typedef std::vector<std::vector<unsigned int>> border_swap_information;

/**
 * @brief Convenience type definition that describes a tuple containing vectors of all flow velocities 
 *        and density values for a fixed time step.
 */
typedef std::tuple<std::vector<velocity>, std::vector<double>> sim_data_tuple;

/**
 * @brief Convenience type definition that describes a tuple containing vectors of all flow velocities 
 *        and density values for a fixed time step. Alternative to sim_data_tuple as the tuple data structure 
 *        sometimes invokes problems.
 */
typedef std::vector<std::vector<velocity>, std::vector<double>> sim_data_vector;

/**
 * @brief This type stands for an access function. Node values can be stored in different layout and 
 *        via this function, the corresponding access scheme can be specified.
 */
typedef std::function<unsigned int(unsigned int, unsigned int)> access_function;

/* Inlet and outlet behavior */
#define INLET_VELOCITY 1.0/100
#define INLET_DENSITY 2
#define OUTLET_DENSITY 2 

/** Mapping of directions as proposed by Mattila to the corresponding velocity vectors */
extern std::map<unsigned int, velocity> velocity_vectors;

/** Mapping of directions as proposed by Mattila to the weights of the corresponding distribution function */
extern std::map<unsigned int,double> weights;

/** Vector containing the distribution values that actually change within a streaming step */
extern std::vector<unsigned int> streaming_directions;

/**
 * @brief Returns the inverse direction of that specified.
 */
inline unsigned int invert_direction(unsigned int dir)
{
    return 8 - dir;
}

/**
 * @brief Returns the Maxwell-Boltzmann-Distribution for all directions in the order proposed by Mattila et al.
 * 
 * @param u two-dimensional velocity vector
 * @param rho density
 * @return the probability of there being a particle with velocity v_direction 
 */
std::vector<double> maxwell_boltzmann_distribution(velocity &u, double rho);

#endif