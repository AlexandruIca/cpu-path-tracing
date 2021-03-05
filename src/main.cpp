#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <fmt/format.h>
#include <taskflow/taskflow.hpp>

#include "random_state.hpp"
#include "ray.hpp"
#include "reflection.hpp"
#include "sphere.hpp"
#include "utils.hpp"
#include "vec.hpp"

constexpr double epsilon = 1e-4;
constexpr double pi = 3.14159265358979323846;

using vec3 = pt::vec3;
using ray = pt::ray;
using reflection_type = pt::reflection_type;
using sphere_t = pt::sphere;

std::array<sphere_t, 10> spheres = { {
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

///
/// \returns The closest intersection point in the whole scene if anything is found, 0 otherwise.
///
[[nodiscard]] auto intersect(const ray& r, double& t, std::size_t& id) noexcept -> double
{
    constexpr double inf = 1e20;
    t = inf;

    for(std::size_t i = 0; i < spheres.size(); i++) {
        if(double const d = spheres.at(i).intersect(r); d > 0 && d < t) {
            t = d;
            id = i;
        }
    }

    return static_cast<double>(t < inf);
}

[[nodiscard]] auto diffuse_ray(vec3 const& hit_point, [[maybe_unused]] vec3 const& normal, pt::rand_state& rng) -> ray
{
    double const phi = 2 * pi * rng.generate();
    double const random_angle = rng.generate(); // 1 - cos^2 theta
    double const sin_theta = std::sqrt(random_angle);
    double const cos_theta = std::sqrt(1.0 - random_angle);

    vec3 const w = normal;
    vec3 const u = (std::abs(w.x) > 0.1 ? vec3{ 0, 1, 0 } : vec3{ 1, 0, 0 }).cross(w).norm();
    vec3 const v = w.cross(u);
    vec3 const new_direction = (u * std::cos(phi) * sin_theta + v * std::sin(phi) * sin_theta + w * cos_theta).norm();

    return ray{ hit_point, new_direction };
}

[[nodiscard]] auto specular_ray(ray const& original, vec3 const& hit_point, vec3 const& outward_normal) -> ray
{
    return ray{ hit_point, original.direction - outward_normal * 2.0 * outward_normal.dot(original.direction) };
}

[[nodiscard]] auto
dielectric_ray(vec3 const& hit_point, vec3 const& uv, vec3 const& normal, double const etai_over_etat) -> ray
{
    auto const cos_theta = std::min((uv * -1.0).dot(normal), 1.0);
    vec3 const r_out_perp = (uv + normal * cos_theta) * etai_over_etat;
    vec3 const r_out_parallel = normal * (-std::sqrt(std::abs(1.0 - r_out_perp.dot(r_out_perp))));

    return ray{ hit_point, r_out_perp + r_out_parallel };
}

[[nodiscard]] auto radiance(const ray& r, int const depth, pt::rand_state& rng) -> vec3
{
    constexpr int russian_roulette_threshold = 4;
    double closest_distance = 0.0;
    std::size_t object_index = 0;

    if(intersect(r, closest_distance, object_index) <= epsilon) {
        return vec3{ 0.0, 0.0, 0.0 };
    }

    const sphere_t& obj = spheres.at(object_index);
    vec3 const hit_point = r.at(closest_distance);
    vec3 const outward_normal = (hit_point - obj.position).norm();
    bool const front_facing = outward_normal.dot(r.direction) < 0;
    // front facing normal:
    vec3 const normal = front_facing ? outward_normal : outward_normal * -1;
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
        return obj.emission + color.blend(radiance(diffuse_ray(hit_point, normal, rng), depth + 1, rng));
    }
    case reflection_type::specular: {
        return obj.emission + color.blend(radiance(specular_ray(r, hit_point, outward_normal), depth + 1, rng));
    }
    case reflection_type::dielectric: {
        constexpr double refraction_index = 2.0;
        double const refraction_ratio = front_facing ? (1.0 / refraction_index) : refraction_index;
        auto ray_in = r;
        auto const unit_direction = ray_in.direction.norm();

        double const cos_theta = std::min((unit_direction * -1.0).dot(normal), 1.0);
        double const sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);

        bool cannot_refract = refraction_ratio * sin_theta > 1.0;
        ray refl{};

        auto reflectance = [](double const cosine, double const ref_idx) noexcept -> double {
            constexpr int offset = 5;
            auto r0 = (1.0 - ref_idx) / (1.0 + ref_idx);
            r0 *= r0;
            return r0 + (1.0 - r0) * std::pow(1.0 - cosine, offset);
        };

        if(cannot_refract || reflectance(cos_theta, refraction_ratio) > rng.generate()) {
            refl = specular_ray(r, hit_point, outward_normal);
        }
        else {
            refl = dielectric_ray(hit_point, unit_direction, normal, refraction_ratio);
        }

        return obj.emission + color.blend(radiance(refl, depth + 1, rng));
    }
    }

    return vec3{ 0, 0, 0 };
}

auto main(int argc, char* argv[]) -> int
{
    std::vector<std::string> args{ argv + 1, argv + argc };
    int constexpr w = 1024;
    int constexpr h = 768;
    int const samps = argc == 2 ? std::stoi(args[0]) / 4 : 1;

    ray cam{ vec3{ 50, 51, 295.6 }, vec3{ 0, -0.042612, -1 }.norm() };
    vec3 const cx = vec3{ w * 0.5135 / h, 0, 0 };
    vec3 const cy = cx.cross(cam.direction).norm() * 0.5135;
    std::vector<vec3> c{};
    c.reserve(w * h);

    tf::Executor executor{};
    tf::Taskflow taskflow{};

    for(int y = 0; y < h; y++) {
        taskflow.emplace([&c, y, samps, cam, cx, cy] {
            std::cerr << fmt::format("\rRendering ({} spp) {:>5.2}%", samps * 4, 100.0 * y / (h - 1));

            auto const seed = static_cast<unsigned short>(y * y * y);
            auto rng = pt::rand_state::default_with_seed(seed);

            for(unsigned short x = 0; x < w; x++) {
                auto const i = static_cast<std::size_t>((h - y - 1) * w) + x;

                // 2x2 subpixel rows
                for(int sy = 0; sy < 2; sy++) {
                    // 2x2 subpixel cols
                    for(int sx = 0; sx < 2; sx++) {
                        vec3 r{ 0, 0, 0 };
                        for(int s = 0; s < samps; s++) {
                            double const r1 = 2.0 * rng.generate();
                            double const dx = r1 < 1.0 ? std::sqrt(r1) - 1.0 : 1.0 - std::sqrt(2.0 - r1);
                            double const r2 = 2.0 * rng.generate();
                            double const dy = r2 < 1.0 ? std::sqrt(r2) - 1.0 : 1.0 - std::sqrt(2.0 - r2);

                            vec3 d = cx * (((sx + 0.5 + dx) / 2 + x) / w - 0.5) +
                                     cy * (((sy + 0.5 + dy) / 2 + y) / h - 0.5) + cam.direction;

                            r = r + radiance(ray{ cam.origin + d * 140, d.norm() }, 0, rng) * (1.0 / samps);
                        }

                        c[i] = c[i] + vec3{ pt::clamp(r.x), pt::clamp(r.y), pt::clamp(r.z) } * 0.25;
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
        g << fmt::format("{} {} {} ", pt::color_to_int(c[i].x), pt::color_to_int(c[i].y), pt::color_to_int(c[i].z));
    }
}
