#ifndef UTILS_HPP
#define UTILS_HPP
#include <array>
#include <vector>
#include <set>
#include <iostream>
#include "defines.hpp"

inline unsigned int matrix_access(unsigned int row, 
    unsigned int column, 
    unsigned int column_count)
{
    return row * column_count + column;
}

template <typename T>
void print_vector(std::vector<T> &vector)
{
    int current_value_count = 0;
    std::cout << "[";
    for(auto i = 0; i < vector.size(); ++i)
    {
        if(current_value_count == HORIZONTAL_NODES)
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

inline void print_phase_vector(std::vector<bool> &vector)
{
    int current_value_count = 0;
    for(auto i = 0; i < vector.size(); ++i)
    {
        if(current_value_count == HORIZONTAL_NODES)
        {
            current_value_count = 0;
            std::cout << std::endl;
        }
        if(vector[i]) std::cout << "#";
        else std::cout << "~"; 
        std::cout << " ";
        current_value_count++;
    }
    std::cout << std::endl;
} 

template <typename T>
void print_vector(std::vector<T> &vector, int row_length)
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

inline void print_velocity_vector(std::vector<velocity> &vector)
{
    int current_value_count = 0;
    std::cout << "[";
    for(auto i = 0; i < vector.size(); ++i)
    {
        if(current_value_count == HORIZONTAL_NODES - 2)
        {
            current_value_count = 0;
            std::cout << std::endl;
        }
        std::cout << "("<< vector[i][0] << ", " << vector[i][1] << ")";
        std::cout << "\t  ";
        current_value_count++;
    }
    std::cout << "]";
    std::cout << std::endl;
} 


template <typename T>
void print_set(std::set<T> &set)
{
    std::cout << "(";
    for(auto element : set)
    {
        std::cout << element;
        std::cout << ", ";
    }
    std::cout << ")";
    std::cout << std::endl;
} 


template <unsigned long int d>
inline void swap(std::array<double, d>& array, unsigned int a, unsigned int b)
{
    double temp = array[a];
    array[a] = array[b];
    array[b] = temp;
}

namespace math_utils
{
    template <unsigned long int d>
    double dot(std::array<double, d> x, std::array<double, d> y)
    {
        double result = 0;
        for(auto i = 0; i < d; ++i) result += x[i] * y[i];
        return result;
    }

    template<unsigned int d>
    double contraction(std::array<double, d> x, std::array<double, d> y)
    {
        double result = 0;
        for(auto i = 0; i < d; ++i) result += x[i] * y[i];
        return result;
    }

    template<unsigned int d>
    std::array<double, d*d> outer(std::array<double, d> x, std::array<double, d> y)
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