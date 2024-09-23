#pragma once

#include "MathFunctions.hpp"
#include "gcem.hpp"
#include <array>
#include <cstddef>

namespace Windows
{

template <std::size_t length>
constexpr auto Kaiser(double const beta) -> std::array<double, length>
{
    std::array<double, length> window{};

    double const denominator = MathFunctions::ModifiedBesselFirstKindZerothOrder(beta);
    double const tmp0 = (static_cast<double>(length) - 1.0) / 2.0;

    for (std::size_t n = 0; n < length; ++n)
    {
        double const tmp1 = (static_cast<double>(n) - tmp0) / tmp0;
        double const tmp2 = beta * gcem::sqrt(1 - tmp1 * tmp1);

        double const numerator = MathFunctions::ModifiedBesselFirstKindZerothOrder(tmp2);

        window[n] = numerator / denominator;
    }

    return window;
}

// #include <numbers>
// constexpr auto Beta(double const attenuation) -> double
// {
//     if (attenuation > 50.0)
//     {
//         return 0.1102 * (attenuation - 8.7);
//     }
//     if (attenuation > 21.0)
//     {
//         return 0.5842 * gcem::pow((attenuation - 21), 0.4) + 0.07886 * (attenuation - 21);
//     }
//     return 0.0;
// }

// constexpr auto FirstNull(size_t const windowLength, double const beta) -> double
// {
//     double const alpha = beta / std::numbers::pi;
//     return gcem::sqrt(1 + alpha * alpha) / static_cast<double>(windowLength);
// }

/*
Compute the attenuation of a Kaiser FIR filter.

    Given the number of taps `N` and the transition width `width`, compute the
    attenuation `a` in dB, given by Kaiser's formula:

        a = 2.285 * (N - 1) * pi * width + 7.95

    # Kaiser's formula (as given in Oppenheim and Schafer) is for the filter
    # order, so we have to add 1 to get the number of taps.
    numtaps = (A - 7.95) / 2.285 / (np.pi * width) + 1



"""Compute the Kaiser parameter `beta`, given the attenuation `a`.

    Parameters
    ----------
    a : float
        The desired attenuation in the stopband and maximum ripple in
        the passband, in dB.  This should be a *positive* number.

    Returns
    -------
    beta : float
        The `beta` parameter to be used in the formula for a Kaiser window.

    References
    ----------
    Oppenheim, Schafer, "Discrete-Time Signal Processing", p.475-476.

    Examples
    --------
    Suppose we want to design a lowpass filter, with 65 dB attenuation
    in the stop band.  The Kaiser window parameter to be used in the
    window method is computed by ``kaiser_beta(65)``:

    >>> from scipy.signal import kaiser_beta
    >>> kaiser_beta(65)
    6.20426

    """
    if a > 50:
        beta = 0.1102 * (a - 8.7)
    elif a > 21:
        beta = 0.5842 * (a - 21) ** 0.4 + 0.07886 * (a - 21)
    else:
        beta = 0.0
    return beta
*/
} // namespace Windows
