#ifndef PT_VEC_HPP
#define PT_VEC_HPP
#pragma once

namespace pt {

struct vec3
{
    double x{ 0.0 };
    double y{ 0.0 };
    double z{ 0.0 };

    vec3() noexcept = delete;
    vec3(vec3 const&) noexcept = default;
    vec3(vec3&&) noexcept = default;
    ~vec3() noexcept = default;

    vec3(double x_, double y_, double z_) noexcept;

    auto operator=(vec3 const&) noexcept -> vec3& = default;
    auto operator=(vec3&&) noexcept -> vec3& = default;

    [[nodiscard]] auto operator+(const vec3& b) const noexcept -> vec3;
    [[nodiscard]] auto operator-(const vec3& b) const noexcept -> vec3;
    [[nodiscard]] auto operator*(double b) const noexcept -> vec3;
    [[nodiscard]] auto blend(const vec3& b) const noexcept -> vec3;
    [[nodiscard]] auto norm() noexcept -> vec3&;
    [[nodiscard]] auto dot(const vec3& b) const noexcept -> double;
    [[nodiscard]] auto cross(vec3 const& b) const noexcept -> vec3;
    [[nodiscard]] auto operator[](int index) const noexcept -> double const&;
};

} // namespace pt

#endif // !PT_VEC_HPP
