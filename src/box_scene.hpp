#ifndef PT_SMALLPT_SCENE_HPP
#define PT_SMALLPT_SCENE_HPP
#pragma once

#include "reflection.hpp"
#include "scene.hpp"
#include "sphere.hpp"
#include "vec.hpp"

#include <array>

namespace pt {

[[nodiscard]] inline auto box_scene(int const w, int const h) -> scene
{
    constexpr double big_radius = 1E6;
    constexpr double offset = 0.4;
    constexpr double y = 0.0;
    constexpr double z = -1.0;

    scene scn{};
    // radius, position, emission, color, material
    scn.spheres.push_back(pt::sphere{ big_radius,
                                      pt::vec3{ -big_radius - offset, y, z },
                                      pt::vec3{ 0.0, 0.0, 0.0 },
                                      pt::vec3{ 0.9, 0.1, 0.2 },
                                      pt::reflection_type::diffuse }); // Left
    scn.spheres.push_back(pt::sphere{ big_radius,
                                      pt::vec3{ big_radius + offset, y, z },
                                      pt::vec3{ 0.0, 0.0, 0.0 },
                                      pt::vec3{ 0.3, 0.1, 0.9 },
                                      pt::reflection_type::diffuse }); // Right
    scn.spheres.push_back(pt::sphere{ big_radius,
                                      pt::vec3{ 0.0, 0.0, z - big_radius },
                                      pt::vec3{ 0.0, 0.0, 0.0 },
                                      pt::vec3{ 0.1, 0.7, 0.2 },
                                      pt::reflection_type::diffuse }); // Back
    scn.spheres.push_back(pt::sphere{ big_radius,
                                      pt::vec3{ 0.0, big_radius + offset, z },
                                      pt::vec3{ 0.0, 0.0, 0.0 },
                                      pt::vec3{ 0.3, 0.7, 0.2 },
                                      pt::reflection_type::diffuse }); // Top
    scn.spheres.push_back(pt::sphere{ big_radius,
                                      pt::vec3{ 0.0, -big_radius - offset, z },
                                      pt::vec3{ 0.0, 0.0, 0.0 },
                                      pt::vec3{ 0.9, 0.9, 0.9 },
                                      pt::reflection_type::diffuse }); // Bottom
    scn.spheres.push_back(pt::sphere{ offset / 2.0,
                                      pt::vec3{ 0.0, 0.0 + offset / 4.0, z - offset / 2.5 },
                                      pt::vec3{ 9.0, 9.0, 9.0 },
                                      pt::vec3{ 1.8, 1.8, 1.8 },
                                      pt::reflection_type::diffuse }); // Light
    scn.spheres.push_back(pt::sphere{ offset / 2.0,
                                      pt::vec3{ offset / 2.0, -offset / 2.0, z + offset * 1.5 },
                                      pt::vec3{ 0.0, 0.0, 0.0 },
                                      pt::vec3{ 1.0, 1.0, 1.0 },
                                      pt::reflection_type::specular }); // Mirror
    scn.spheres.push_back(pt::sphere{ offset / 2.0,
                                      pt::vec3{ -offset / 2.0, -offset / 2.0, z + offset * 1.5 },
                                      pt::vec3{ 0.0, 0.0, 0.0 },
                                      pt::vec3{ 1.0, 1.0, 1.0 },
                                      pt::reflection_type::dielectric }); // Dielectric

    scn.camera_parameters.position = pt::vec3{ 0.0, 0.0, 2.0 };
    scn.camera_parameters.direction = pt::vec3{ 0.0, 0.0, z + offset * 1.5 };
    scn.camera_parameters.aspect_ratio = (w * 1.0) / (h * 1.0);
    scn.camera_parameters.vertical_fov_radians = 0.5;
    scn.camera_parameters.aperture = 0.2;
    scn.camera_parameters.focus_distance = (scn.camera_parameters.position - scn.camera_parameters.direction).length();

    return scn;
}

} // namespace pt

#endif // !PT_SMALLPT_SCENE_HPP
