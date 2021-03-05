#ifndef PT_SMALLPT_SCENE_HPP
#define PT_SMALLPT_SCENE_HPP
#pragma once

#include "reflection.hpp"
#include "sphere.hpp"
#include "vec.hpp"

#include <array>

namespace pt {

using sphere_t = sphere;

inline std::array<sphere_t, 10> spheres = { {
    // Scene: radius, position, emission, color, material
    sphere_t{ 1e5,
              vec3{ 1e5 + 1, 40.8, 81.6 },
              vec3{ 0.0, 0.0, 0.0 },
              vec3{ 0.75, 0.25, 0.25 },
              reflection_type::diffuse }, // Left
    sphere_t{ 1e5,
              vec3{ -1e5 + 99, 40.8, 81.6 },
              vec3{ 0.0, 0.0, 0.0 },
              vec3{ 0.25, 0.25, 0.75 },
              reflection_type::diffuse }, // Right
    sphere_t{
        1e5, vec3{ 50, 40.8, 1e5 }, vec3{ 0.0, 0.0, 0.0 }, vec3{ 0.75, 0.75, 0.75 }, reflection_type::diffuse }, // Back
    sphere_t{ 1e5,
              vec3{ 50, 40.8, -1e5 + 170 },
              vec3{ 0.0, 0.0, 0.0 },
              vec3{ 0.0, 0.0, 0.0 },
              reflection_type::diffuse }, // Front
    sphere_t{ 1e5,
              vec3{ 50, 1e5, 81.6 },
              vec3{ 0.0, 0.0, 0.0 },
              vec3{ 0.75, 0.75, 0.75 },
              reflection_type::diffuse }, // Bottom
    sphere_t{ 1e5,
              vec3{ 50, -1e5 + 81.6, 81.6 },
              vec3{ 0.0, 0.0, 0.0 },
              vec3{ 0.25, 0.75, 0.15 },
              reflection_type::diffuse }, // Top
    sphere_t{ 16.5,
              vec3{ 27, 16.5, 47 },
              vec3{ 0.0, 0.0, 0.0 },
              vec3{ 1, 1, 1 } * 0.999,
              reflection_type::specular }, // Mirror
    sphere_t{ 16.5,
              vec3{ 65, 16.5, 37 },
              vec3{ 0.0, 0.0, 0.0 },
              vec3{ 0.6, 0.1, 0.6 },
              reflection_type::specular }, // Mirror Purple
    sphere_t{
        16.5, vec3{ 45, 46.5, 50 }, vec3{ 22, 22, 22 }, vec3{ 0.0, 0.0, 0.0 }, reflection_type::diffuse }, // Light up
    sphere_t{ 16.5,
              vec3{ 73, 16.5, 78 },
              vec3{ 0.0, 0.0, 0.0 },
              vec3{ 1, 1, 1 } * 0.999,
              reflection_type::dielectric }, // Glass
} };

} // namespace pt

#endif // !PT_SMALLPT_SCENE_HPP
