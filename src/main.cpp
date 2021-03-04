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

constexpr double epsilon = 1e-4;
constexpr double pi = 3.14159265358979323846;

struct vec3
{
    double x{ 0.0 };
    double y{ 0.0 };
    double z{ 0.0 };

    vec3() noexcept = delete;
    vec3(vec3 const&) noexcept = default;
    vec3(vec3&&) noexcept = default;
    ~vec3() noexcept = default;

    vec3(double x_, double y_, double z_) noexcept
        : x{ x_ }
        , y{ y_ }
        , z{ z_ }
    {
    }

    auto operator=(vec3 const&) noexcept -> vec3& = default;
    auto operator=(vec3&&) noexcept -> vec3& = default;

    [[nodiscard]] auto operator+(const vec3& b) const noexcept -> vec3
    {
        return vec3{ x + b.x, y + b.y, z + b.z };
    }

    [[nodiscard]] auto operator-(const vec3& b) const noexcept -> vec3
    {
        return vec3{ x - b.x, y - b.y, z - b.z };
    }

    [[nodiscard]] auto operator*(double b) const noexcept -> vec3
    {
        return vec3{ x * b, y * b, z * b };
    }

    [[nodiscard]] auto mult(const vec3& b) const noexcept -> vec3
    {
        return vec3{ x * b.x, y * b.y, z * b.z };
    }

    [[nodiscard]] auto norm() noexcept -> vec3&
    {
        return *this = *this * (1 / sqrt(x * x + y * y + z * z));
    }

    [[nodiscard]] auto dot(const vec3& b) const noexcept -> double
    {
        return x * b.x + y * b.y + z * b.z;
    }

    [[nodiscard]] auto cross(vec3 const& b) const noexcept -> vec3
    {
        return vec3{ y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x };
    }
};

struct Ray
{
    vec3 o{ 0, 0, 0 };
    vec3 d{ 0, 0, 0 };
};

enum Refl_t
{
    DIFF,
    SPEC,
    REFR
};

struct Sphere
{
    double rad{ 0.0 };
    vec3 p{ 0, 0, 0 };
    vec3 e{ 0, 0, 0 };
    vec3 c{ 0, 0, 0 };
    Refl_t refl{ Refl_t::DIFF };

    ///
    /// \returns 0 if no intersection was found, something greater than 0 otherwise
    ///
    [[nodiscard]] auto intersect(const Ray& r) const noexcept -> double
    {
        vec3 op = p - r.o;
        double t = 0.0;
        double const b = op.dot(r.d);
        double det = b * b - op.dot(op) + rad * rad;

        if(det < 0) {
            return 0;
        }

        det = sqrt(det);
        t = b - det;

        if(t > epsilon) {
            return t;
        }

        t = b + det;

        if(t > epsilon) {
            return t;
        }

        return 0.0;
    }
};

std::array<Sphere, 10> spheres = { {
    // Scene: radius, position, emission, color, material
    Sphere{ 1e5, vec3{ 1e5 + 1, 40.8, 81.6 }, vec3{ 0.0, 0.0, 0.0 }, vec3{ 0.75, 0.25, 0.25 }, DIFF },   // Left
    Sphere{ 1e5, vec3{ -1e5 + 99, 40.8, 81.6 }, vec3{ 0.0, 0.0, 0.0 }, vec3{ 0.25, 0.25, 0.75 }, DIFF }, // Right
    Sphere{ 1e5, vec3{ 50, 40.8, 1e5 }, vec3{ 0.0, 0.0, 0.0 }, vec3{ 0.75, 0.75, 0.75 }, DIFF },         // Back
    Sphere{ 1e5, vec3{ 50, 40.8, -1e5 + 170 }, vec3{ 0.0, 0.0, 0.0 }, vec3{ 0.0, 0.0, 0.0 }, DIFF },     // Front
    Sphere{ 1e5, vec3{ 50, 1e5, 81.6 }, vec3{ 0.0, 0.0, 0.0 }, vec3{ 0.75, 0.75, 0.75 }, DIFF },         // Bottom
    Sphere{ 1e5, vec3{ 50, -1e5 + 81.6, 81.6 }, vec3{ 0.0, 0.0, 0.0 }, vec3{ 0.25, 0.75, 0.15 }, DIFF }, // Top
    Sphere{ 16.5, vec3{ 27, 16.5, 47 }, vec3{ 0.0, 0.0, 0.0 }, vec3{ 1, 1, 1 } * 0.999, SPEC },          // Mirror
    Sphere{ 16.5, vec3{ 65, 16.5, 37 }, vec3{ 0.0, 0.0, 0.0 }, vec3{ 0.6, 0.1, 0.6 }, SPEC },   // Mirror Purple
    Sphere{ 16.5, vec3{ 45, 46.5, 50 }, vec3{ 22, 22, 22 }, vec3{ 0.0, 0.0, 0.0 }, DIFF },      // Light up
    Sphere{ 16.5, vec3{ 73, 16.5, 78 }, vec3{ 0.0, 0.0, 0.0 }, vec3{ 1, 1, 1 } * 0.999, REFR }, // Glass
    // Sphere(600, Vec(50, 681.6, 0.27, 81.6), Vec(6, 6, 6), Vec(0.2, 0.2, 0.5), DIFF) // Light
} };

///
/// \brief Clamps x to (0.0, 1.0).
///
[[nodiscard]] auto clamp(double const x) noexcept -> double
{
    return std::clamp(x, 0.0, 1.0);
}

[[nodiscard]] auto toInt(double const x) noexcept -> int
{
    // gamma 2.2 correction
    double const corrected = std::pow(clamp(x), 1.0 / 2.2);
    return static_cast<int>(std::round(corrected * 255.0));
}

///
/// \returns The closest intersection point in the whole scene if anything is found, 0 otherwise.
///
[[nodiscard]] auto intersect(const Ray& r, double& t, std::size_t& id) noexcept -> double
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

[[nodiscard]] auto diffuse_ray(vec3 const& hit_point, vec3 const& normal, pt::rand_state& rng) -> Ray
{
    double const r1 = 2 * pi * rng.generate(); // phi
    double const r2 = rng.generate();          // 1 - cos^2 theta
    double const r2s = sqrt(r2);               // sin_theta

    vec3 const w = normal;
    vec3 const u = (fabs(w.x) > 0.1 ? vec3{ 0, 1, 0 } : vec3{ 1, 0, 0 }).cross(w).norm();
    vec3 const v = w.cross(u);
    vec3 const d = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).norm();

    return Ray{ hit_point, d };
    // double const r1 = erand48(Xi);
    // double const r2 = erand48(Xi);
    // double const sin_phi = sqrt(r1); // r1 == 1 - cos^2(phi)
    // double const cos_phi = sqrt(std::max(0.0, 1 - r1));
    // double const theta = 2 * 3.14159265 * r2;

    // auto const new_direction = Vec{ sin_phi * cos(theta), sin_phi * sin(theta), cos_phi };

    // return obj.e + f.mult(radiance(Ray(x, new_direction), depth, Xi));
}

[[nodiscard]] auto specular_ray(Ray const& original, vec3 const& hit_point, vec3 const& outward_normal) -> Ray
{
    return Ray{ hit_point, original.d - outward_normal * 2.0 * outward_normal.dot(original.d) };
}

[[nodiscard]] auto
dielectric_ray(vec3 const& hit_point, vec3 const& uv, vec3 const& normal, double const etai_over_etat) -> Ray
{
    auto const cos_theta = std::min((uv * -1.0).dot(normal), 1.0);
    vec3 const r_out_perp = (uv + normal * cos_theta) * etai_over_etat;
    vec3 const r_out_parallel = normal * (-std::sqrt(std::abs(1.0 - r_out_perp.dot(r_out_perp))));

    return Ray{ hit_point, r_out_perp + r_out_parallel };
}

[[nodiscard]] auto radiance(const Ray& r, int const depth, pt::rand_state& rng) -> vec3
{
    double t = 0.0;
    std::size_t id = 0;

    if(intersect(r, t, id) <= epsilon) {
        return vec3{ 0.0, 0.0, 0.0 };
    }

    const Sphere& obj = spheres.at(id);
    vec3 const x = r.o + r.d * t;
    vec3 const n = (x - obj.p).norm();
    vec3 const nl = n.dot(r.d) < 0 ? n : n * -1;
    vec3 f = obj.c;

    double const p = std::max({ f.x, f.y, f.z });

    if(depth > 4) {
        // Russina Roulette
        if(rng.generate() < p) {
            f = f * (1.0 / p);
        }
        else {
            return obj.e;
        }
    }

    switch(obj.refl) {
    case DIFF: {
        return obj.e + f.mult(radiance(diffuse_ray(x, nl, rng), depth + 1, rng));
    }
    case SPEC: {
        return obj.e + f.mult(radiance(specular_ray(r, x, n), depth + 1, rng));
    }
    case REFR: {
        constexpr double refraction_index = 2.0;
        double const refraction_ratio = (n.dot(nl) > 0) ? (1.0 / refraction_index) : refraction_index;
        auto ray_in = r;
        auto const unit_direction = ray_in.d.norm();

        double const cos_theta = std::min((unit_direction * -1.0).dot(nl), 1.0);
        double const sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);

        bool cannot_refract = refraction_ratio * sin_theta > 1.0;
        Ray refl{};

        auto reflectance = [](double const cosine, double const ref_idx) noexcept -> double {
            constexpr int offset = 5;
            auto r0 = (1.0 - ref_idx) / (1.0 + ref_idx);
            r0 *= r0;
            return r0 + (1.0 - r0) * std::pow(1.0 - cosine, offset);
        };

        if(cannot_refract || reflectance(cos_theta, refraction_ratio) > rng.generate()) {
            refl = specular_ray(r, x, n);
        }
        else {
            refl = dielectric_ray(x, unit_direction, nl, refraction_ratio);
        }

        return obj.e + f.mult(radiance(refl, depth + 1, rng));
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

    Ray cam{ vec3{ 50, 51, 295.6 }, vec3{ 0, -0.042612, -1 }.norm() };
    vec3 const cx = vec3{ w * 0.5135 / h, 0, 0 };
    vec3 const cy = cx.cross(cam.d).norm() * 0.5135;
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
                            double const dx = r1 < 1.0 ? sqrt(r1) - 1.0 : 1.0 - sqrt(2.0 - r1);
                            double const r2 = 2.0 * rng.generate();
                            double const dy = r2 < 1.0 ? sqrt(r2) - 1.0 : 1.0 - sqrt(2.0 - r2);

                            vec3 d = cx * (((sx + 0.5 + dx) / 2 + x) / w - 0.5) +
                                     cy * (((sy + 0.5 + dy) / 2 + y) / h - 0.5) + cam.d;

                            r = r + radiance(Ray{ cam.o + d * 140, d.norm() }, 0, rng) * (1.0 / samps);
                        }

                        c[i] = c[i] + vec3{ clamp(r.x), clamp(r.y), clamp(r.z) } * 0.25;
                    }
                }
            }
        });
    }

    executor.run(taskflow).wait();

    std::ofstream g{ "image.ppm" };
    g << fmt::format("P3\n{} {}\n{}\n", w, h, 255);

    for(std::size_t i = 0; i < w * h; ++i) {
        g << fmt::format("{} {} {} ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z));
    }
}
