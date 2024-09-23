#pragma once

#include "MathFunctions.hpp"
#include <array>
#include <cstddef>

namespace Filters
{
template <std::size_t length>
constexpr auto LowPass(
    std::array<double, length> const window,
    double const normalisedCutoffFrequency) -> std::array<double, length>
{
    std::array<double, length> filter{};

    auto const offset = (static_cast<double>(length) - 1.0) / 2.0;
    auto const scale = 2.0 * normalisedCutoffFrequency;
    for (std::size_t i = 0; i < length; ++i)
    {
        auto const t = static_cast<double>(i) - offset;
        filter[i] = window[i] * scale * MathFunctions::Sinc(scale * t);
    }

    return filter;
}
} // namespace Filters