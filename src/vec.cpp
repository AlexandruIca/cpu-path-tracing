#include "vec.hpp"

#include <cmath>

namespace pt {

vec3::vec3(double const x_, double const y_, double const z_) noexcept
    : x{ x_ }
    , y{ y_ }
    , z{ z_ }
{
}

auto vec3::operator+(vec3 const& b) const noexcept -> vec3
{
    return vec3{ x + b.x, y + b.y, z + b.z };
}

auto vec3::operator-(vec3 const& b) const noexcept -> vec3
{
    return vec3{ x - b.x, y - b.y, z - b.z };
}

auto vec3::operator*(double const b) const noexcept -> vec3
{
    return vec3{ x * b, y * b, z * b };
}

auto vec3::blend(vec3 const& b) const noexcept -> vec3
{
    return vec3{ x * b.x, y * b.y, z * b.z };
}

auto vec3::norm() noexcept -> vec3&
{
    return *this = *this * (1 / std::sqrt(x * x + y * y + z * z));
}

auto vec3::dot(vec3 const& b) const noexcept -> double
{
    return x * b.x + y * b.y + z * b.z;
}

auto vec3::cross(vec3 const& b) const noexcept -> vec3
{
    return vec3{ y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x };
}

} // namespace pt
