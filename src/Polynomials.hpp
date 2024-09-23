#pragma once

#include <array>
#include <cstddef>

namespace Polynomials::Chebyshev
{
template <typename T, std::size_t length>
constexpr auto Evaluate(std::array<T, length> const& coefficients, T const x) -> T
{
    if constexpr (length == 0)
    {
        return static_cast<T>(0);
    }
    else if constexpr (length == 1)
    {
        return coefficients[0];
    }
    else if constexpr (length == 2)
    {
        return coefficients[0] + coefficients[1] * x;
    }
    else
    {
        T const x2 = x * static_cast<T>(2);

        T t0 = coefficients[length - 2];
        T t1 = coefficients[length - 1];

        for (size_t i = 3; i < length + 1; ++i)
        {
            T const tmp = t0;
            t0 = coefficients[length - i] - t1;
            t1 = x2 * t1 + tmp;
        }

        return t0 + t1 * x;
    }
}
} // namespace Polynomials::Chebyshev

// #include <span>
// template <typename T>
// constexpr auto Evaluate(std::span<T> const& coefficients, T const x) -> T
// {
//     if (coefficients.size() == 0)
//     {
//         return static_cast<T>(0);
//     }

//     if (coefficients.size() == 1)
//     {
//         return coefficients[0];
//     }

//     if (coefficients.size() == 2)
//     {
//         return coefficients[0] + coefficients[1] * x;
//     }

//     std::size_t const length = coefficients.size();
//     T const x2 = x * static_cast<T>(2);

//     T t0 = coefficients[length - 2];
//     T t1 = coefficients[length - 1];

//     for (size_t i = 3; i < length + 1; ++i)
//     {
//         T const tmp = t0;
//         t0 = coefficients[length - i] - t1;
//         t1 = x2 * t1 + tmp;
//     }

//     return t0 + t1 * x;
// }

// Evaluate polynomials that are of the first order Chebyshev basis. Draws from scipy's chebval
// implementation.
// template <typename T>
// constexpr auto Evaluate(std::array<T, 1> const& coefficients, T const /*x*/) -> T
// {
//     return coefficients[0];
// }

// template <typename T>
// constexpr auto Evaluate(std::array<T, 2> const& coefficients, T const x) -> T
// {
//     return coefficients[0] + coefficients[1] * x;
// }

// template <typename T, std::size_t length>
// constexpr auto Evaluate(std::array<T, length> const& coefficients, T const x) -> T
// {
//     static_assert(
//         length > 2,
//         "Failed to evaluated Chebyshev Polynomial because there were insufficient
//         coefficients.");

//     T const x2 = x * static_cast<T>(2);

//     T t0 = coefficients[length - 2];
//     T t1 = coefficients[length - 1];

//     for (size_t i = 3; i < length + 1; ++i)
//     {
//         T const tmp = t0;
//         t0 = coefficients[length - i] - t1;
//         t1 = x2 * t1 + tmp;
//     }

//     return t0 + t1 * x;
// }

// template <typename T>
// constexpr auto Evaluate(std::span<T> const& coefficients, T const x) -> T
// {
//     if (coefficients.size() == 0)
//     {
//         return 0;
//     }

//     if (coefficients.size() == 1)
//     {
//         return coefficients[0];
//     }

//     if (coefficients.size() == 2)
//     {
//         return coefficients[0] + coefficients[1] * x;
//     }

//     std::size_t const length = coefficients.size();
//     T const x2 = x * static_cast<T>(2);

//     T t0 = coefficients[length - 2];
//     T t1 = coefficients[length - 1];

//     for (size_t i = 3; i < length + 1; ++i)
//     {
//         T const tmp = t0;
//         t0 = coefficients[length - i] - t1;
//         t1 = x2 * t1 + tmp;
//     }

//     return t0 + t1 * x;
// }
// } // namespace Polynomials