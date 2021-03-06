#include <algorithm>
#include <array>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <fmt/format.h>
#include <taskflow/taskflow.hpp>

#include "constants.hpp"
#include "hit_record.hpp"
#include "random_state.hpp"
#include "ray.hpp"
#include "reflection.hpp"
#include "smallpt_scene.hpp"
#include "sphere.hpp"
#include "utils.hpp"
#include "vec.hpp"

using vec3 = pt::vec3;
using ray = pt::ray;
using reflection_type = pt::reflection_type;

///
/// \returns The closest intersection point in the whole scene if anything is found, 0 otherwise.
///
[[nodiscard]] auto intersect(const ray& r, double& t, std::size_t& id) noexcept -> bool
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

[[nodiscard]] auto diffuse_ray(pt::hit_record const& record, pt::rand_state& rng) -> ray
{
    double const phi = 2 * pt::pi * rng.generate();
    double const random_angle = rng.generate(); // 1 - cos^2 theta
    double const sin_theta = std::sqrt(random_angle);
    double const cos_theta = std::sqrt(1.0 - random_angle);

    vec3 const w = record.normal;
    vec3 const u = (std::abs(w.x) > 0.1 ? vec3{ 0, 1, 0 } : vec3{ 1, 0, 0 }).cross(w).norm();
    vec3 const v = w.cross(u);
    vec3 const new_direction = (u * std::cos(phi) * sin_theta + v * std::sin(phi) * sin_theta + w * cos_theta).norm();

    return ray{ record.hit_point, new_direction };
}

[[nodiscard]] auto specular_ray(pt::hit_record const& record, pt::rand_state& rng) -> ray
{
    double constexpr fuzziness = 0.0;
    auto const reflected = record.original_ray.direction -
                           record.outward_normal * 2.0 * record.outward_normal.dot(record.original_ray.direction);
    auto const factor = rng.generate() * fuzziness;
    return ray{ record.hit_point, reflected + vec3{ factor, factor, factor } };
}

[[nodiscard]] auto dielectric_ray(pt::hit_record const& record, pt::rand_state& rng) -> ray
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

    vec3 const r_out_perp = (unit_direction + record.normal * cos_theta) * refraction_ratio;
    vec3 const r_out_parallel = record.normal * (-std::sqrt(std::abs(1.0 - r_out_perp.dot(r_out_perp))));

    return ray{ record.hit_point, r_out_perp + r_out_parallel };
}

[[nodiscard]] auto radiance(ray const& r, int const depth, pt::rand_state& rng) -> vec3
{
    constexpr int russian_roulette_threshold = 4;
    double closest_distance = 0.0;
    std::size_t object_index = 0;

    if(!intersect(r, closest_distance, object_index)) {
        return vec3{ 0.0, 0.0, 0.0 };
    }

    auto const& obj = pt::spheres.at(object_index);
    auto const record = get_hit_record_at(obj, r, closest_distance);
    vec3 color = obj.color;

    double const probability = std::max({ color.x, color.y, color.z });

    if(depth > russian_roulette_threshold) {
        if(rng.generate() < probability) {
            color = color * (1.0 / probability);
        }
        else {
            return obj.emission;
        }
    }

    switch(obj.reflection) {
    case reflection_type::diffuse: {
        return obj.emission + color.blend(radiance(diffuse_ray(record, rng), depth + 1, rng));
    }
    case reflection_type::specular: {
        return obj.emission + color.blend(radiance(specular_ray(record, rng), depth + 1, rng));
    }
    case reflection_type::dielectric: {
        return obj.emission + color.blend(radiance(dielectric_ray(record, rng), depth + 1, rng));
    }
    }

    return vec3{ 0, 0, 0 };
}

auto main(int argc, char* argv[]) -> int
{
    // how many subpixels per row/column?
    constexpr int num_subpixels = 2;
    std::vector<std::string> args{ argv + 1, argv + argc };
    int constexpr w = 1024;
    int constexpr h = 768;
    int const samps = argc == 2 ? std::stoi(args[0]) / num_subpixels : 1;

    ray cam{ vec3{ 50, 51, 295.6 }, vec3{ 0, -0.042612, -1 }.norm() };

    constexpr double aspect_ratio = (w * 1.0) / (h * 1.0);
    constexpr double vertical_fov = 0.502643; // approx. 28 degrees
    double const fov_scale = 2.0 * std::tan(0.5 * vertical_fov);

    vec3 const cam_x_axis = vec3{ aspect_ratio, 0, 0 } * fov_scale;
    // cam_y_axis is orghogonal to both x and forward axis
    vec3 const cam_y_axis = cam_x_axis.cross(cam.direction).norm() * fov_scale;
    std::vector<vec3> image{};
    image.resize(w * h, vec3{ 0, 0, 0 });

    tf::Executor executor{};
    tf::Taskflow taskflow{};

    for(int y = 0; y < h; y++) {
        taskflow.emplace([&image, y, samps, cam, cam_x_axis, cam_y_axis] {
            std::cerr << fmt::format(
                "\rRendering ({} spp) {:>5.2}%", samps * num_subpixels * num_subpixels, 100.0 * y / (h - 1));

            auto const seed = static_cast<unsigned short>(y * y * y);
            auto rng = pt::rand_state::default_with_seed(seed);

            for(unsigned short x = 0; x < w; x++) {
                auto const i = static_cast<std::size_t>((h - y - 1) * w) + x;

                for(int sy = 0; sy < num_subpixels; sy++) {
                    for(int sx = 0; sx < num_subpixels; sx++) {
                        vec3 r{ 0, 0, 0 };

                        for(int s = 0; s < samps; s++) {
                            constexpr double subpixel_length = 1.0 / num_subpixels;
                            double const x_in_subpixel = (x + sx * subpixel_length + rng.generate() * subpixel_length);
                            double const y_in_subpixel = (y + sy * subpixel_length + rng.generate() * subpixel_length);

                            // At what point on the x/y axis are we on? (between -0.5 and 0.5)
                            vec3 const offset_x = cam_x_axis * (x_in_subpixel / w - 0.5);
                            vec3 const offset_y = cam_y_axis * (y_in_subpixel / h - 0.5);
                            vec3 new_direction = cam.direction + offset_x + offset_y;

                            ray const new_ray{ cam.origin + new_direction * 140.0, new_direction.norm() };
                            vec3 const contributor = radiance(new_ray, 0, rng);
                            r = r + contributor * (1.0 / samps);
                        }

                        vec3 const subpixel_color = vec3{ pt::clamp(r.x), pt::clamp(r.y), pt::clamp(r.z) };
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
        g << fmt::format(
            "{} {} {} ", pt::color_to_int(image[i].x), pt::color_to_int(image[i].y), pt::color_to_int(image[i].z));
    }
}
