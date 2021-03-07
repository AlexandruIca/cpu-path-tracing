#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <fmt/format.h>
#include <taskflow/taskflow.hpp>

#include "camera.hpp"
#include "constants.hpp"
#include "hit_record.hpp"
#include "random_state.hpp"
#include "ray.hpp"
#include "reflection.hpp"
#include "sphere.hpp"
#include "utils.hpp"
#include "vec.hpp"

#include "simple_scene.hpp"

///
/// \returns The closest intersection point in the whole scene if anything is found, 0 otherwise.
///
[[nodiscard]] auto intersect(const pt::ray& r, double& t, std::size_t& id) noexcept -> bool
{
    t = pt::inf;

    for(std::size_t i = 0; i < pt::spheres.size(); i++) {
        if(double const d = pt::spheres.at(i).intersect(r); d > 0 && d < t) {
            t = d;
            id = i;
        }
    }

    return t < pt::inf;
}

[[nodiscard]] auto diffuse_ray(pt::hit_record const& record, pt::rand_state& rng) -> pt::ray
{
    double const phi = 2 * pt::pi * rng.generate();
    double const random_angle = rng.generate(); // 1 - cos^2 theta
    double const sin_theta = std::sqrt(random_angle);
    double const cos_theta = std::sqrt(1.0 - random_angle);

    pt::vec3 const w = record.normal;
    pt::vec3 const u = (std::abs(w.x) > 0.1 ? pt::vec3{ 0, 1, 0 } : pt::vec3{ 1, 0, 0 }).cross(w).norm();
    pt::vec3 const v = w.cross(u);
    pt::vec3 const new_direction =
        (u * std::cos(phi) * sin_theta + v * std::sin(phi) * sin_theta + w * cos_theta).norm();

    return pt::ray{ record.hit_point, new_direction };
}

[[nodiscard]] auto specular_ray(pt::hit_record const& record, pt::rand_state& rng) -> pt::ray
{
    double constexpr fuzziness = 0.0;
    auto const reflected = record.original_ray.direction -
                           record.outward_normal * 2.0 * record.outward_normal.dot(record.original_ray.direction);
    auto const factor = rng.generate() * fuzziness;
    return pt::ray{ record.hit_point, reflected + pt::vec3{ factor, factor, factor } };
}

[[nodiscard]] auto dielectric_ray(pt::hit_record const& record, pt::rand_state& rng) -> pt::ray
{
    constexpr double refraction_index = 2.0;
    double const refraction_ratio = record.front_facing ? (1.0 / refraction_index) : refraction_index;

    auto r = record.original_ray;
    auto const unit_direction = r.direction.norm();

    double const cos_theta = std::min((unit_direction * -1.0).dot(record.normal), 1.0);
    double const sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);

    bool cannot_refract = refraction_ratio * sin_theta > 1.0;

    auto reflectance = [](double const cosine, double const ref_idx) noexcept -> double {
        constexpr int offset = 5;
        auto r0 = (1.0 - ref_idx) / (1.0 + ref_idx);
        r0 *= r0;
        return r0 + (1.0 - r0) * std::pow(1.0 - cosine, offset);
    };

    if(cannot_refract || reflectance(cos_theta, refraction_ratio) > rng.generate()) {
        return specular_ray(record, rng);
    }

    pt::vec3 const r_out_perp = (unit_direction + record.normal * cos_theta) * refraction_ratio;
    pt::vec3 const r_out_parallel = record.normal * (-std::sqrt(std::abs(1.0 - r_out_perp.dot(r_out_perp))));

    return pt::ray{ record.hit_point, r_out_perp + r_out_parallel };
}

///
/// emisson + reflectance * radiance
/// Iteratively:
/// emisson0 + reflectance0 * emisson1 + reflectance0 * reflectance1 * emisson2 + ...
///
[[nodiscard]] auto radiance(pt::ray const& ray, pt::rand_state& rng) -> pt::vec3
{
    constexpr int russian_roulette_threshold = 4;
    pt::vec3 accumulated_emission{ 0, 0, 0 };
    pt::vec3 accumulated_reflectance{ 1, 1, 1 };
    pt::ray r = ray;

    for(int depth = 0; depth < pt::depth_limit; ++depth) {
        double closest_distance = 0.0;
        std::size_t object_index = 0;

        if(!intersect(r, closest_distance, object_index)) {
            return accumulated_emission;
        }

        auto const& obj = pt::spheres.at(object_index);
        auto const record = get_hit_record_at(obj, r, closest_distance);
        pt::vec3 color = obj.color;

        accumulated_emission = accumulated_emission + accumulated_reflectance.blend(obj.emission);

        double const probability = std::max({ color.x, color.y, color.z });

        if(depth > russian_roulette_threshold) {
            if(rng.generate() < probability) {
                color = color * (1.0 / probability);
            }
            else {
                return accumulated_emission;
            }
        }

        accumulated_reflectance = accumulated_reflectance.blend(color);

        switch(obj.reflection) {
        case pt::reflection_type::diffuse: {
            r = diffuse_ray(record, rng);
            break;
        }
        case pt::reflection_type::specular: {
            r = specular_ray(record, rng);
            break;
        }
        case pt::reflection_type::dielectric: {
            r = dielectric_ray(record, rng);
            break;
        }
        }
    }

    return accumulated_emission;
}

auto main(int argc, char* argv[]) -> int
{
    // how many subpixels per row/column?
    constexpr int num_subpixels = 2;
    std::vector<std::string> args{ argv + 1, argv + argc };
    int constexpr w = 1024;
    int constexpr h = 768;
    [[maybe_unused]] int const samps = argc == 2 ? std::stoi(args[0]) / num_subpixels : 1;

    pt::camera_config cfg{};
    cfg.position = pt::vec3{ -2.0, 2.0, 1.0 };
    cfg.direction = pt::vec3{ 0.0, 0.0, -1.0 };
    cfg.aspect_ratio = (w * 1.0) / (h * 1.0);
    cfg.vertical_fov_radians = 1.2;
    cfg.aperture = 0.5;
    cfg.focus_distance = (cfg.position - cfg.direction).length();

    auto const cam = pt::camera::with_config(cfg);
    std::vector<pt::vec3> image{};
    image.resize(w * h, pt::vec3{ 0, 0, 0 });

    tf::Executor executor{};
    tf::Taskflow taskflow{};

    for(int y = 0; y < h; y++) {
        taskflow.emplace([&image, y, samps, &cam] {
            std::cerr << fmt::format(
                "\rRendering ({} spp) {:>5.2}%", samps * num_subpixels * num_subpixels, 100.0 * y / (h - 1));

            auto const seed = static_cast<unsigned short>(y * y * y);
            auto rng = pt::rand_state::default_with_seed(seed);

            for(unsigned short x = 0; x < w; x++) {
                auto const i = static_cast<std::size_t>((h - y - 1) * w) + x;

                for(int sy = 0; sy < num_subpixels; sy++) {
                    for(int sx = 0; sx < num_subpixels; sx++) {
                        pt::vec3 r{ 0, 0, 0 };

                        for(int s = 0; s < samps; s++) {
                            // At what point on the x/y axis are we on? (between 0.0 and 1.0)
                            constexpr double subpixel_length = 1.0 / num_subpixels;
                            double const x_in_subpixel =
                                (x + sx * subpixel_length + rng.generate_between(0.0, subpixel_length));
                            double const y_in_subpixel =
                                (y + sy * subpixel_length + rng.generate_between(0.0, subpixel_length));

                            pt::ray const new_ray{ cam.get_ray(x_in_subpixel / w, y_in_subpixel / h, rng) };
                            pt::vec3 const contributor = radiance(new_ray, rng);
                            r = r + contributor * (1.0 / samps);
                        }

                        pt::vec3 const subpixel_color = pt::vec3{ pt::clamp(r.x), pt::clamp(r.y), pt::clamp(r.z) };
                        image[i] = image[i] + subpixel_color * (1.0 / (num_subpixels * num_subpixels));
                    }
                }
            }
        });
    }

    executor.run(taskflow).wait();

    std::cerr << std::endl;

    std::ofstream g{ "image.ppm" };
    g << fmt::format("P3\n{} {}\n{}\n", w, h, 255);

    for(std::size_t i = 0; i < w * h; ++i) {
        g << fmt::format("{} ", pt::color_to_int(image[i].x));
        g << fmt::format("{} ", pt::color_to_int(image[i].y));
        g << fmt::format("{} ", pt::color_to_int(image[i].z));
    }
}
