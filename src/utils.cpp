#include "utils.hpp"

#include <algorithm>
#include <cmath>

auto pt::clamp(double const x) noexcept -> double
{
    return std::clamp(x, 0.0, 1.0);
}

auto pt::color_to_int(double const x) noexcept -> int
{
    // gamma 2.2 correction
    double const corrected = std::pow(clamp(x), 1.0 / 2.2);
    return static_cast<int>(std::round(corrected * 255.0));
}
