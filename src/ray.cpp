#include "ray.hpp"

auto pt::ray::at(double const t) const noexcept -> pt::vec3
{
    return origin + direction * t;
}
