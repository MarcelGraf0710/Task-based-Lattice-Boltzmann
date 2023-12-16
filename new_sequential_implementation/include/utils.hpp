#ifndef UTILS_HPP
#define UTILS_HPP
#include <array>

inline unsigned int matrix_access(unsigned int row, 
    unsigned int column, 
    unsigned int column_count)
{
    return row * column_count + column;
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