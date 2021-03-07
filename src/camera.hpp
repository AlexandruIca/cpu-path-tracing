#ifndef PT_CAMERA_HPP
#define PT_CAMERA_HPP
#pragma once

#include "random_state.hpp"
#include "ray.hpp"
#include "vec.hpp"

namespace pt {

struct camera_config
{
    vec3 position{ 0, 0, 0 };
    vec3 direction{ 0, 0, 0 };
    vec3 up{ 0, 1, 0 };
    double aspect_ratio{ 16.0 / 9.0 };
    double vertical_fov_radians{ 0.785398163 }; // default value is approximately 45 degrees
    double focal_length{ 1.0 };
    double aperture{ 0.0 };
    double focus_distance{ 0.0 };
};

struct camera
{
    vec3 position{ 0, 0, 0 };
    vec3 lower_left_corner{ 0, 0, 0 };
    vec3 cam_x_axis{ 0, 0, 0 };
    vec3 cam_y_axis{ 0, 0, 0 };
    vec3 u{ 0, 0, 0 };
    vec3 v{ 0, 0, 0 };
    vec3 w{ 0, 0, 0 };
    double lens_radius{ 0.0 };

    [[nodiscard]] static auto with_config(camera_config const& cfg) noexcept -> camera;
    [[nodiscard]] static auto random_in_unit_disk(pt::rand_state& rng) -> vec3;
    ///
    /// \param u The x coordinate between 0.0 and 1.0.
    /// \param v The y coordinate between 0.0 and 1.0.
    ///
    /// Basically divide x by width and y by height before passing them here.
    ///
    [[nodiscard]] auto get_ray(double s, double t, pt::rand_state& rng) const noexcept -> ray;
};

} // namespace pt

#endif // !PT_CAMERA_HPP
