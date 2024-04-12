#ifndef DEFINES_HPP
#define DEFINES_HPP

/// Include directives ///

#include <vector>
#include <list>
#include <complex>
#include <map>
#include <functional>
#include <valarray>
#include <numeric>

/// General definitions (NOTE: this may be outsourced to a .csv file in the future to increase modifiability.) ///

#define DIMENSION_COUNT 2
#define DIRECTION_COUNT 9
#define VERTICAL_NODES 26 // 15, readable up to 50 in console, beyond not recommended // 14 as for parallel
#define HORIZONTAL_NODES 7 // comfortable console limit: 9 
#define TOTAL_NODE_COUNT VERTICAL_NODES * HORIZONTAL_NODES
#define BOLTZMANN_CONSTANT 1.380649e-23
#define RELAXATION_TIME 1.4
#define TIME_STEPS 50
#define SUBDOMAIN_HEIGHT 8 // 4 as for parallel
#define SUBDOMAIN_COUNT (VERTICAL_NODES + 1) / (1 + SUBDOMAIN_HEIGHT)
#define BUFFER_COUNT SUBDOMAIN_COUNT - 1

/// Inlet and outlet behavior ///

#define INLET_VELOCITY velocity{0.1,0.0}
#define OUTLET_VELOCITY velocity{0.00,0.0}
#define INLET_DENSITY 1
#define OUTLET_DENSITY 1 

/// Convenience and readability type definitions ///

/**
 * @brief Representation of a velocity vector
 */
typedef std::array<double, DIMENSION_COUNT> velocity; 

/**
 * @brief Convenience type definition that represents a vector from which the boundary treatment 
 *        of all nodes can be retrieved. Information is stored in the following way:
 *        Each entry of the outer vector represents a border node, called BORDER_NODE for explanation.
 *        Every BORDER_NODE is a vector with the following information:
 *        - 0th entry: The index of BORDER_NODE
 *        - Further entries: The directions pointing to non-inout ghost nodes (including solid nodes within the domain).
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
 * @brief This type stands for an access function. Node values can be stored in different layout and 
 *        via this function, the corresponding access scheme can be specified.
 */
typedef std::function<unsigned int(unsigned int, unsigned int)> access_function;

/// Global constants ///

/** Mapping of directions as proposed by Mattila to the corresponding velocity vectors */
extern const std::map<unsigned int, velocity> VELOCITY_VECTORS;

/** Mapping of directions as proposed by Mattila to the WEIGHTS of the corresponding distribution function */
extern const std::map<unsigned int,double> WEIGHTS;

/** Vector containing the distribution values that actually change within a streaming step */
extern const std::vector<unsigned int> STREAMING_DIRECTIONS;

/** Vector containing all directions */
extern const std::vector<unsigned int> ALL_DIRECTIONS;

/// Methods ///

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
std::vector<double> maxwell_boltzmann_distribution(const velocity &u, const double rho);

#endif