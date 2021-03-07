#include "camera.hpp"

[[nodiscard]] auto pt::camera::with_config(camera_config const& cfg) noexcept -> camera
{
    double const viewport_height = 2.0 * std::tan(0.5 * cfg.vertical_fov_radians);
    double const viewport_width = cfg.aspect_ratio * viewport_height;

    auto const w = (cfg.position - cfg.direction).norm();
    auto const u = cfg.up.cross(w).norm();
    auto const v = w.cross(u);

    vec3 const cam_x_axis = u * viewport_width * cfg.focus_distance;
    vec3 const cam_y_axis = v * viewport_height * cfg.focus_distance;
    vec3 const lower_left_corner = cfg.position - cam_x_axis * 0.5 - cam_y_axis * 0.5 - w * cfg.focus_distance;

    return camera{ cfg.position, lower_left_corner, cam_x_axis, cam_y_axis, u, v, w, cfg.aperture / 2.0 };
}

[[nodiscard]] auto pt::camera::random_in_unit_disk(pt::rand_state& rng) -> vec3
{
    for(;;) {
        vec3 const point{ rng.generate_between(-1.0, 1.0), rng.generate_between(-1.0, 1.0), 0.0 };

        if(point.dot(point) >= 1.0) {
            continue;
        }

        return point;
    }
}

[[nodiscard]] auto pt::camera::get_ray(double const s, double const t, pt::rand_state& rng) const noexcept -> ray
{
    auto const rd = random_in_unit_disk(rng) * lens_radius;
    vec3 const offset{ rd * s + rd * t };
    vec3 direction{ lower_left_corner + cam_x_axis * s + cam_y_axis * t - position - offset };
    return ray{ position + offset, direction };
}
