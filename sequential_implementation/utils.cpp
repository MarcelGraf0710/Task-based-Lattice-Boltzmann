#include "utils.hpp"

template <unsigned long int d>
double math_utils::dot(std::array<double, d> x, std::array<double, d> y)
{
    double result = 0;
    for(auto i = 0; i < d; ++i) result += x[i] * y[i];
    return result;
}

template<unsigned int d>
double math_utils::contraction(std::array<double, d> x, std::array<double, d> y)
{
    double result = 0;
    for(auto i = 0; i < d; ++i) result += x[i] * y[i];
    return result;
}

template<unsigned int d>
std::array<double, d*d> math_utils::outer(std::array<double, d> x, std::array<double, d> y)
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

