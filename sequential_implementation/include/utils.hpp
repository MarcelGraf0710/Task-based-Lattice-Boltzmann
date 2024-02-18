#ifndef UTILS_HPP
#define UTILS_HPP
#include <array>
#include <vector>
#include <set>
#include <iostream>
#include "defines.hpp"
#include "access.hpp"
#include "boundaries.hpp"
#include <iomanip>

inline unsigned int matrix_access
(
    unsigned int row, 
    unsigned int column, 
    unsigned int column_count)
{
    return row * column_count + column;
}

/**
 * @brief This namespace contains multiple utility functions that allow the user
 *        to print a data structure to the console.
 */
namespace to_console
{
    /**
     * @brief Allows to print out a vector representing a data layout whose column count is HORIZONTAL_NODES.
     *        Notice that the vector is assumed to represent a matrix.
     * 
     * @tparam T the type of the objects the specified vector holds (must be numeric)
     * @param vector the vector that is to be printed in the console
     */
    template <typename T>
    void print_vector(const std::vector<T> &vector)
    {
        for(auto y = VERTICAL_NODES - 1; y >= 0; --y)
        {
            for(auto x = 0; x < HORIZONTAL_NODES; ++x)
            {
                if(x == 0 && y == 0) std::cout << "\033[31m";
                else if(x == (HORIZONTAL_NODES - 1) && y == (VERTICAL_NODES -1)) std::cout << "\033[34m";
                std::cout << vector[matrix_access(y,x, HORIZONTAL_NODES)];
                std::cout << "\t\033[0m";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    } 

    /**
     * @brief Adapted version of print_vector that supports a custom column count.
     *        In order to enforce that the vector is printed as a row vector, specify a column count
     *        that exceeds the length of the vector, e.g. std::numeric_limits<int>::max_value()
     * 
     * @tparam T the type of the objects the specified vector holds (must be numeric)
     * @param vector the vector that is to be printed in the console
     * @param row_length a line break will be set after this many entries
     */
    template <typename T>
    void print_vector(const std::vector<T> &vector, unsigned int row_length)
    {
        int current_value_count = 0;
        std::cout << "[";
        for(auto i = 0; i < vector.size(); ++i)
        {
            if(current_value_count == row_length)
            {
                current_value_count = 0;
                std::cout << std::endl;
            }
            std::cout << vector[i];
            std::cout << "\t";
            current_value_count++;
        }
        std::cout << "]";
        std::cout << std::endl;
    }

    /**
     * @brief Adapted version of print_vector that prints out the phase vector of a lattice.
     *        If a node is solid (i.e. the entry is true), it is represented by #.
     *        If a node is fluid (i.e. the entry is false), it is represented by ~.
     * 
     * @param vector the phase vector
     */
    inline void print_phase_vector(const std::vector<bool> &vector)
    {
        for(auto y = VERTICAL_NODES - 1; y >= 0; --y)
        {
            for(auto x = 0; x < HORIZONTAL_NODES; ++x)
            {
                if(vector[matrix_access(y,x, HORIZONTAL_NODES)]) std::cout << "\033[32m#\033[0m";
                else std::cout << "\033[34m~\033[0m"; 
                std::cout << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    } 

    /**
     * @brief Prints all velocity values in the lattice to the console.
     *        All values are printed in order, i.e. the origin is located at the lower left corner of the output.
     * 
     * @param vector a vector containing all velocity values
     */
    inline void print_velocity_vector(const std::vector<velocity> &vector)
    {
        for(auto y = VERTICAL_NODES - 1; y >= 0; --y)
        {
            for(auto x = 0; x < HORIZONTAL_NODES; ++x)
            {
                if(x == 0 && y == 0) std::cout << "\033[31m";
                else if(x == (HORIZONTAL_NODES - 1) && y == (VERTICAL_NODES -1)) std::cout << "\033[34m";
                std::cout << "("<< vector[matrix_access(y,x, HORIZONTAL_NODES)][0] << ", " << vector[matrix_access(y,x, HORIZONTAL_NODES)][1] << ")";
                std::cout << "\t  \033[0m";
                std::cout << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    } 

    /**
     * @brief Prints the contents of a set to the console.
     * 
     * @tparam T the type of the objects the set holds (must be compatible with the console)
     * @param set the set whose contents are to be printed
     */
    template <typename T>
    void print_set(const std::set<T> &set)
    {
        std::cout << "(";
        for(const auto element : set)
        {
            std::cout << element;
            std::cout << ", ";
        }
        std::cout << ")";
        std::cout << std::endl;
    }

    /**
     * @brief Prints all distribution values in to the console.
     *        They are displayed in the original order, i.e. the origin is located at the lower left corner of the printed distribution chart.
     * 
     * @param distribution_values a vector containing the distribution values of all nodes 
     * @param access_function the function used to access the distribution values
     */
    inline void print_distribution_values
    (
        const std::vector<double> &distribution_values, 
        const access_function access_function
    )
    {
        std::vector<std::vector<unsigned int>> print_dirs = {{6,7,8}, {3,4,5}, {0,1,2}};
        unsigned int current_node_index = 0;
        unsigned int previous_direction = 0;
        std::vector<double> current_values(9,0);
        std::cout << std::setprecision(3) << std::fixed;

        for(auto y = VERTICAL_NODES - 1; y >= 0; --y)
        {
            for(auto i = 0; i < 3; ++i)
            {
                auto current_row = print_dirs[i];
                for(auto x = 0; x < HORIZONTAL_NODES; ++x)
                {
                    if(x == 0 && y == 0) std::cout << "\033[31m";
                    else if(x == (HORIZONTAL_NODES - 1) && y == (VERTICAL_NODES -1)) std::cout << "\033[34m";
                    current_node_index = access::get_node_index(x, y);
                    current_values = access::get_distribution_values_of(distribution_values, current_node_index, access_function);

                    for(auto j = 0; j < 3; ++j)
                    {
                        auto direction = current_row[j];
                        std::cout << current_values[direction] << "  ";
                    }
                    std::cout << "\t\033[0m";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
            std::cout << std::endl;
        }
    } 
}

/**
 * @brief This namespace contains some utility functions for vectors.
 */
namespace vec_utils
{
    /**
     * @brief Swaps the entries of the specified vector.
     */
    inline void swap(std::vector<double> &vector, unsigned int a, unsigned int b)
    {
        double temp = vector[a];
        vector[a] = vector[b];
        vector[b] = temp;
    }
}

/**
 * @brief This namespace contains some functions that are convenient for mathematical computations.
 */
namespace math_utils
{
    /**
     * @brief Calculates the dot product of the specified arrays
     */
    template <unsigned long int d>
    double dot(const std::array<double, d> &x, const std::array<double, d> &y)
    {
        double result = 0;
        for(auto i = 0; i < d; ++i) result += x[i] * y[i];
        return result;
    }

    /**
     * @brief Calculates the contraction of two matrices stored as arrays in row-major order.
     */
    template<unsigned int d>
    double contraction(const std::array<double, d> &x, const std::array<double, d> &y)
    {
        double result = 0;
        for(auto i = 0; i < d; ++i) result += x[i] * y[i];
        return result;
    }

    /**
     * @brief Calculates the outer product of the specified two arrays.
     * @return the resulting matrix stored as an array in row-major order.
     */
    template<unsigned int d>
    std::array<double, d*d> outer(const std::array<double, d> &x, const std::array<double, d> &y)
    {
        std::array<double, d*d> result;
        for(auto i = 0; i < d; ++i)
        {
            for(auto j = 0; j < d; ++j)
            {
                result[matrix_access(i,j,d)] = x[i] * y[j];
            }
        }
        return result;
    }
}

#endif