#ifndef PT_SMALLPT_SCENE_HPP
#define PT_SMALLPT_SCENE_HPP
#pragma once

#include "reflection.hpp"
#include "sphere.hpp"
#include "vec.hpp"

#include <array>

namespace pt {

inline std::array<sphere, 5> spheres = { {
    // Scene: radius, position, emission, color, material
    sphere{ 100.0,
            vec3{ 0.0, -100.5, -1.0 },
            vec3{ 0.0, 0.0, 0.0 },
            vec3{ 0.8, 0.8, 0.0 },
            reflection_type::diffuse }, // Ground
    sphere{ 0.5,
            vec3{ 1.0, 0.0, -1.0 },
            vec3{ 0.0, 0.0, 0.0 },
            vec3{ 0.999, 0.999, 0.999 },
            reflection_type::specular }, // Right
    sphere{ 0.5,
            vec3{ -1.0, 0.0, -1.0 },
            vec3{ 0.0, 0.0, 0.0 },
            vec3{ 0.999, 0.999, 0.999 },
            reflection_type::dielectric }, // left
    sphere{ 0.5,
            vec3{ 0.0, 0.0, -1.0 },
            vec3{ 0.1, 0.1, 0.9 },
            vec3{ 0.0, 0.7, 0.1 },
            reflection_type::diffuse }, // Light center
    sphere{ 1.0,
            vec3{ 1.0, 3.1, -1.0 },
            vec3{ 30.0, 30.0, 30.0 },
            vec3{ 0.0, 0.0, 0.0 },
            reflection_type::diffuse }, // Light up
} };

} // namespace pt

#endif // !PT_SMALLPT_SCENE_HPP
