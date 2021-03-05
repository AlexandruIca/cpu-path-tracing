#ifndef PT_RAY_HPP
#define PT_RAY_HPP
#pragma once

#include "vec.hpp"

namespace pt {

struct ray
{
    vec3 origin{ 0, 0, 0 };
    vec3 direction{ 0, 0, 0 };

    [[nodiscard]] auto at(double t) const noexcept -> vec3;
};

} // namespace pt

#endif // !PT_RAY_HPP
