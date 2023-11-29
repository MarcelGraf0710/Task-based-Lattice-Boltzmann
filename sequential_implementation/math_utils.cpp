#include <array>

inline unsigned int matrix_access(unsigned int row, 
    unsigned int column, 
    unsigned int column_count)
{
    return row * column_count + column;
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
        for(i = 0; i < d; ++i)
        {
            for(j = 0; j < d; ++j)
            {
                result[matrix_access(i,j,d)] = x[i] * y[j];
            }
        }
        return result;
    }
}
