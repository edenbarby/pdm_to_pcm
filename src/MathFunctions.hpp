#pragma once

#include "Polynomials.hpp"
#include "gcem.hpp"
#include <array>

namespace MathFunctions
{
constexpr auto ModifiedBesselFirstKindZerothOrder(double const x) -> double
{
    double result = 0.0;
    double const x0 = gcem::abs(x);

    constexpr double firstPolynomialExtent = 8.0;
    if (x0 <= firstPolynomialExtent)
    {
        double const x1 = (x0 / 4.0) - 1.0;
        constexpr std::array<double, 30> coefficients
            = {0.33839763720473803,    -0.3046826723431984,     0.17162090152220877,
               -0.09490109704804764,   0.04930528423967071,     -0.02373741480589947,
               0.010546460394594998,   -0.004324309995050576,   0.0016394756169413357,
               -0.0005763755745385824, 0.00018850288509584165,  -5.754195010082104e-05,
               1.6448448070728896e-05, -4.4167383584587505e-06, 1.1173875391201037e-06,
               -2.670793853940612e-07, 6.046995022541919e-08,   -1.300025009986248e-08,
               2.6598237246823866e-09, -5.189795601635263e-10,  9.675809035373237e-11,
               -1.726826291441556e-11, 2.95505266312964e-12,    -4.856446783111929e-13,
               7.676185498604936e-14,  -1.1685332877993451e-14, 1.715391285555133e-15,
               -2.431279846547955e-16, 3.3307945188222384e-17,  -4.4153416464793395e-18};
        result = gcem::exp(x0) * Polynomials::Chebyshev::Evaluate(coefficients, x1);
    }
    else
    {
        double const x1 = (16.0 / x0) - 1.0;
        constexpr std::array<double, 25> coefficients
            = {0.4022452055070544,      0.0033691164782556943,   6.889758346916825e-05,
               2.8913705208347567e-06,  2.0489185894690638e-07,  2.266668990498178e-08,
               3.3962320257083865e-09,  4.94060238822497e-10,    1.1889147107846439e-11,
               -3.1499165279632416e-11, -1.3215811840447713e-11, -1.7941785315068062e-12,
               7.180124451383666e-13,   3.8527783827421426e-13,  1.54008621752141e-14,
               -4.150569347287222e-14,  -9.554846698828307e-15,  3.8116806693526224e-15,
               1.7725601330565263e-15,  -3.425485619677219e-16,  -2.8276239805165836e-16,
               3.461222867697461e-17,   4.46562142029676e-17,    -4.830504485944182e-18,
               -7.233180487874754e-18};
        result
            = gcem::exp(x0) * Polynomials::Chebyshev::Evaluate(coefficients, x1) / gcem::sqrt(x0);
    }

    return result;
}

constexpr auto Sinc(double const x) -> double
{
    if (x == 0)
    {
        return 1;
    }
    return gcem::sin(x) / x;
}
} // namespace MathFunctions