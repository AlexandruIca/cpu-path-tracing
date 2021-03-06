#ifndef PT_HIT_RECORD_HPP
#define PT_HIT_RECORD_HPP
#pragma once

#include "ray.hpp"
#include "sphere.hpp"
#include "vec.hpp"

namespace pt {

struct hit_record
{
    ray original_ray{};
    vec3 hit_point{ 0, 0, 0 };
    vec3 outward_normal{ 0, 0, 0 };
    vec3 normal{ 0, 0, 0 };
    bool front_facing{ false };
};

[[nodiscard]] auto get_hit_record_at(sphere const& sphere, ray const& r, double t) noexcept -> hit_record;

} // namespace pt

#endif // !PT_HIT_RECORD_HPP
