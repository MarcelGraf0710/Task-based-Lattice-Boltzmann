/* Include directives */
#include <vector>
#include <list>
#include "gauss_hermite_quadrature.cpp"
#include <complex>
#include "math_utils.cpp"

/* General definitions */
#define DIMENSION_COUNT 2
#define DIRECTION_COUNT 9
#define TIME_STEP 1
#define SPACE_STEP 1
#define MOLECULAR_MASS 1
#define VERTICAL_NODES 100
#define HORIZONTAL_NODES 100
#define TOTAL_NODE_COUNT VERTICAL_NODES * HORIZONTAL_NODES
#define SPEED_OF_SOUND 340 
#define TEMPERATURE 20
#define BOLTZMANN_CONSTANT 1.380649e-23
#define RELAXATION_FREQUENCY 1

/* Velocity vectors */
std::map<int, std::array<double,2>> velocity_vectors 
{
    {6, {-1, 1}},  {7, {0, -1}}, {8, {1, 1}},   
    {3, {-1, 0}},  {4, {0, 0}},  {5, {1, 0}},   
    {0, {-1, -1}}, {1, {0, -1}}, {2, {1, -1}}   
};

/**
 * @brief This namespace contains functions that map input values to array index accesses.
 */
namespace access
{
    /**
     * @brief Returns the array index of the collision layout
     * 
     * @param node the node in the simulation domain
     * @param direction the direction of the velocity vector
     * @return the index of the array storing the distribution values 
     */
    inline unsigned int collision(unsigned int node, unsigned int direction)
    {
        return DIRECTION_COUNT * node + direction;        
    }

    /**
     * @brief Returns the array index of the stream layout
     * 
     * @param node the node in the simulation domain
     * @param direction the direction of the velocity vector
     * @return the index of the array storing the distribution values  
     */
    inline unsigned int stream(unsigned int node, unsigned int direction)
    {
        return DIRECTION_COUNT * direction + node;
    }

    /**
     * @brief Returns the array index of the bundle layout
     * 
     * @param node the node in the simulation domain
     * @param direction the direction of the velocity vector
     * @return the index of the array storing the distribution values  
     */
    inline unsigned int bundle(unsigned int node, unsigned int direction)
    {
        return (direction / 3) * TOTAL_NODE_COUNT + (direction % 3); 
    }

    /**
     * @brief Returns the index the desired node has within the array that stores it. 
     *        The origin lies at the lower left corner and enumeration is row-major.
     * 
     * @param x x coordinate
     * @param y y coordinate
     * @return the index of the desired note
     */
    inline unsigned int node(unsigned int x, unsigned int y)
    {
        return x + y * HORIZONTAL_NODES;
    }

}

/**
 * @brief This namespace contains all function representations of boundary conditions used in the lattice-Boltzmann model.
 */
namespace boundary_conditions
{

}

/**
 * @brief This namespace contains all functions used in the lattice-Boltzmann model.
 */
namespace distribution_functions
{

    /**
     * @brief The Maxwell-Boltzmann-Distribution marks the equilibrium distribution of particles.
     * 
     * @param u two-dimensional velocity vector
     * @param rho density
     * @param direction direction according to the scheme proposed by Mattila et al.
     * @return the probability of there being a particle with velocity v_direction 
     */
    double maxwell_boltzmann_distribution (std::array<double, 2> u, double rho, unsigned int direction)
    {
        /*
        std::array<int, DIMENSION_COUNT> c_i = velocity_vectors[direction];
        double quadratic = (1 / pow(SPEED_OF_SOUND, 2)) * u[0] * c_i[0];
        double power_four_term = (1 / 2 * pow(SPEED_OF_SOUND, 4)) *  
                                    (
                                        pow(u[0], 2) * (pow(c_i[0], 2) - 1) 
                                        + 2 * u[0] * u[1] * c_i[0] * c_i[1]
                                        + pow(u[1], 2) * (pow(c_i[1], 2))
                                    );
        return weights[direction] * rho * quadratic * power_four_term;
        */
       return weights[direction] * rho * (1 + 3 * math_utils::dot(velocity_vectors[direction], u) 
        + 9.0/2 * pow(math_utils::dot(velocity_vectors[direction], u), 2) - 3.0/2 * math_utils::dot(velocity_vectors[direction], u));
    }
}

namespace semi_direct_access
{
    /**
     * @brief Returns a vector containing all fluid segments in the specified domain. The order is as follows:
     *        All even numbers mark the beginning of a fluid node.
     *        All odd numbers state how many fluid tiles there are in a row, i.e. how many consecutive fluid tiles there are.
     *        For example, the sequence "3, 10" means that the nodes 3...12 are fluid nodes.
     * 
     * @param phase_space a vector containing the phase information for each node where "true" means solid
     * @return a vector containing the fluid segments in the explained arrangement
     */
    std::vector<unsigned int> get_fluid_segments(std::vector<bool> phase_space)
    {
        unsigned int index = 0;
        unsigned int consecution = 0;
        std::list<unsigned int> fluid_segments;

        while (index < size(phase_space))
        {
            if (!phase_space[index + consecution]) consecution++; // Hit fluid node
            else if (consecution > 0) // Hit solid node that marks the end of a line of consecutive fluid nodes
            {
                fluid_segments.push_back(index);
                fluid_segments.push_back(consecution);
                index += consecution;
                consecution = 0;
            }
            else index++; // Hit consecutive solid nodes
        }

        std::vector<unsigned int> result{
            std::make_move_iterator(begin(fluid_segments)), 
            std::make_move_iterator(end(fluid_segments))
            };

        return result;
    }
}

