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
#include "scene.hpp"
#include "utils.hpp"
#include "vec.hpp"

//#include "simple_scene.hpp"
//#include "box_scene.hpp"
#include "box_mirror_scene.hpp"

///
/// \returns The closest intersection point in the whole scene if anything is found, 0 otherwise.
///
[[nodiscard]] auto intersect(pt::scene const& scene, const pt::ray& r, double& t, std::size_t& id) noexcept -> bool
{
    t = pt::inf;

    for(std::size_t i = 0; i < scene.spheres.size(); i++) {
        if(double const d = scene.spheres.at(i).intersect(r); d > 0 && d < t) {
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
[[nodiscard]] auto radiance(pt::scene const& scene, pt::ray const& ray, pt::rand_state& rng) -> pt::vec3
{
    constexpr int russian_roulette_threshold = 4;
    pt::vec3 accumulated_emission{ 0.0, 0.0, 0.0 };
    pt::vec3 accumulated_reflectance{ 1, 1, 1 };
    pt::ray r = ray;

    for(int depth = 0; depth < pt::depth_limit; ++depth) {
        double closest_distance = 0.0;
        std::size_t object_index = 0;

        if(!intersect(scene, r, closest_distance, object_index)) {
            auto const unit_direction = r.direction.norm();
            auto const t = 0.5 * (unit_direction.y + 1.0);
            auto const background = pt::vec3{ 1.0, 1.0, 1.0 } * (1.0 - t) + pt::vec3{ 0.5, 0.7, 1.0 } * t;
            return accumulated_emission + accumulated_reflectance.blend(background);
        }

        auto const& obj = scene.spheres.at(object_index);
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

struct render_state
{
    pt::scene const& scene;
    pt::camera const& cam;
    pt::rand_state& rng;
    std::vector<pt::vec3>& image;
    int width = 0;
    int height = 0;
    int num_samples = 0;
    int num_subpixels = 0;
};

///
/// \param x Where are we on the x axis.
/// \param y Where are we on the y axis.
/// \param sx what fraction of pixel on the x axis.
/// \param sy what fraction of pixel on the y axis.
/// \param state See [render_state](\ref render_state).
///
auto render_subpixel(int const x, int const y, int const sx, int const sy, render_state& state) -> void
{
    auto const row = static_cast<std::size_t>((state.height - y - 1) * state.width) + static_cast<std::size_t>(x);
    pt::vec3 r{ 0, 0, 0 };

    for(int s = 0; s < state.num_samples; s++) {
        // At what point on the x/y axis are we on? (between 0.0 and 1.0)
        double const subpixel_length = 1.0 / state.num_subpixels;
        double const x_in_subpixel = (x + sx * subpixel_length + state.rng.generate_between(0.0, subpixel_length));
        double const y_in_subpixel = (y + sy * subpixel_length + state.rng.generate_between(0.0, subpixel_length));

        pt::ray const new_ray = state.cam.get_ray(x_in_subpixel / state.width, y_in_subpixel / state.height, state.rng);
        pt::vec3 const contributor = radiance(state.scene, new_ray, state.rng);
        r = r + contributor * (1.0 / state.num_samples);
    }

    pt::vec3 const subpixel_color = pt::vec3{ pt::clamp(r.x), pt::clamp(r.y), pt::clamp(r.z) };
    state.image[row] = state.image[row] + subpixel_color * (1.0 / (state.num_subpixels * state.num_subpixels));
}

auto main(int argc, char* argv[]) -> int
{
    // how many subpixels per row/column?
    constexpr int num_subpixels = 2;
    std::vector<std::string> args{ argv + 1, argv + argc };
    int constexpr width = 1024;
    int constexpr height = 768;
    int const samps = argc == 2 ? std::stoi(args[0]) / (num_subpixels * num_subpixels) : 1;

    auto const some_scene = pt::box_scene(width, height);
    auto const cam = pt::camera::with_config(some_scene.camera_parameters);
    std::vector<pt::vec3> image{};

    image.resize(width * height, pt::vec3{ 0, 0, 0 });

    tf::Executor executor{};
    tf::Taskflow taskflow{};

    for(int y = 0; y < height; y++) {
        taskflow.emplace([&image, y, samps, &cam, &some_scene] {
            std::cerr << fmt::format(
                "\rRendering ({} spp) {:>5.2}%", samps * num_subpixels * num_subpixels, 100.0 * y / (height - 1));

            auto const seed = static_cast<unsigned short>(y * y * y);
            auto rng = pt::rand_state::default_with_seed(seed);
            render_state state{ some_scene, cam, rng, image, width, height, samps, num_subpixels };

            for(unsigned short x = 0; x < width; x++) {
                for(int sy = 0; sy < num_subpixels; sy++) {
                    for(int sx = 0; sx < num_subpixels; sx++) {
                        render_subpixel(x, y, sx, sy, state);
                    }
                }
            }
        });
    }

    executor.run(taskflow).wait();

    std::cerr << std::endl;

    std::ofstream g{ "image.ppm" };
    g << fmt::format("P3\n{} {}\n{}\n", width, height, 255);

    for(std::size_t i = 0; i < width * height; ++i) {
        g << fmt::format("{} ", pt::color_to_int(image[i].x));
        g << fmt::format("{} ", pt::color_to_int(image[i].y));
        g << fmt::format("{} ", pt::color_to_int(image[i].z));
    }
}
