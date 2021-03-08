#ifndef PT_SMALLPT_SCENE_HPP
#define PT_SMALLPT_SCENE_HPP
#pragma once

#include "reflection.hpp"
#include "scene.hpp"
#include "sphere.hpp"
#include "vec.hpp"

#include <array>

namespace pt {

[[nodiscard]] inline auto simple_scene(int const w, int const h) -> scene
{
    scene scn{};
    // radius, position, emission, color, material
    scn.spheres.push_back(pt::sphere{ 100.0,
                                      pt::vec3{ 0.0, -100.5, -1.0 },
                                      pt::vec3{ 0.0, 0.0, 0.0 },
                                      pt::vec3{ 0.8, 0.8, 0.0 },
                                      pt::reflection_type::diffuse }); // Ground
    scn.spheres.push_back(pt::sphere{ 0.5,
                                      pt::vec3{ 1.0, 0.0, -1.0 },
                                      pt::vec3{ 0.0, 0.0, 0.0 },
                                      pt::vec3{ 0.999, 0.999, 0.999 },
                                      pt::reflection_type::specular }); // Right
    scn.spheres.push_back(pt::sphere{ 0.5,
                                      pt::vec3{ -1.0, 0.0, -1.0 },
                                      pt::vec3{ 0.0, 0.0, 0.0 },
                                      pt::vec3{ 0.999, 0.999, 0.999 },
                                      pt::reflection_type::dielectric }); // left
    scn.spheres.push_back(pt::sphere{ 0.5,
                                      pt::vec3{ 0.0, 0.0, -1.0 },
                                      pt::vec3{ 0.1, 0.1, 0.9 },
                                      pt::vec3{ 0.0, 0.7, 0.1 },
                                      pt::reflection_type::diffuse }); // Light center
    scn.spheres.push_back(pt::sphere{ 1.0,
                                      pt::vec3{ 1.0, 3.1, -1.0 },
                                      pt::vec3{ 30.0, 30.0, 30.0 },
                                      pt::vec3{ 0.0, 0.0, 0.0 },
                                      pt::reflection_type::diffuse }); // Light up

    scn.camera_parameters.position = pt::vec3{ -2.0, 2.0, 1.0 };
    scn.camera_parameters.direction = pt::vec3{ 0.0, 0.0, -1.0 };
    scn.camera_parameters.aspect_ratio = (w * 1.0) / (h * 1.0);
    scn.camera_parameters.vertical_fov_radians = 1.2;
    scn.camera_parameters.aperture = 0.2;
    scn.camera_parameters.focus_distance = (scn.camera_parameters.position - scn.camera_parameters.direction).length();

    return scn;
}

} // namespace pt

#endif // !PT_SMALLPT_SCENE_HPP
