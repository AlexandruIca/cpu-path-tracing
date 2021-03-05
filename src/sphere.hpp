#ifndef PT_SPHERE_HPP
#define PT_SPHERE_HPP
#pragma once

#include "ray.hpp"
#include "reflection.hpp"

namespace pt {

struct sphere
{
    double radius{ 0.0 };
    vec3 position{ 0, 0, 0 };
    vec3 emission{ 0, 0, 0 };
    vec3 color{ 0, 0, 0 };
    reflection_type reflection{ reflection_type::diffuse };

    ///
    /// \returns 0 if no intersection was found, something greater than 0 otherwise
    ///
    [[nodiscard]] auto intersect(ray const& r) const noexcept -> double;
};

} // namespace pt

#endif // !PT_SPHERE_HPP
