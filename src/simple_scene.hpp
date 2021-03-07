#ifndef PT_SMALLPT_SCENE_HPP
#define PT_SMALLPT_SCENE_HPP
#pragma once

#include "reflection.hpp"
#include "sphere.hpp"
#include "vec.hpp"

#include <array>

namespace pt {

using sphere_t = sphere;

inline std::array<sphere_t, 5> spheres = { {
    // Scene: radius, position, emission, color, material
    sphere_t{ 100.0,
              vec3{ 0.0, -100.5, -1.0 },
              vec3{ 0.0, 0.0, 0.0 },
              vec3{ 0.8, 0.8, 0.0 },
              reflection_type::diffuse }, // Ground
    sphere_t{ 0.5,
              vec3{ 1.0, 0.0, -1.0 },
              vec3{ 0.0, 0.0, 0.0 },
              vec3{ 0.999, 0.999, 0.999 },
              reflection_type::specular }, // Right
    sphere_t{ 0.5,
              vec3{ -1.0, 0.0, -1.0 },
              vec3{ 0.0, 0.0, 0.0 },
              vec3{ 0.999, 0.999, 0.999 },
              reflection_type::dielectric }, // left
    sphere_t{ 0.5,
              vec3{ 0.0, 0.0, -1.0 },
              vec3{ 0.1, 0.1, 0.9 },
              vec3{ 0.0, 0.7, 0.1 },
              reflection_type::diffuse }, // Light center
    sphere_t{ 1.0,
              vec3{ 1.0, 3.0, -1.0 },
              vec3{ 30.0, 30.0, 30.0 },
              vec3{ 0.0, 0.0, 0.0 },
              reflection_type::diffuse }, // Light up
} };

} // namespace pt

#endif // !PT_SMALLPT_SCENE_HPP
