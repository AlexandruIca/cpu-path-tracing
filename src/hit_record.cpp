#include "hit_record.hpp"

[[nodiscard]] auto pt::get_hit_record_at(sphere const& sphere, ray const& r, double const t) noexcept -> hit_record
{
    vec3 const hit_point = r.at(t);
    vec3 const outward_normal = (hit_point - sphere.position).norm();
    bool const front_facing = outward_normal.dot(r.direction) < 0;
    // front facing normal:
    vec3 const normal = front_facing ? outward_normal : outward_normal * -1;

    return { r, hit_point, outward_normal, normal, front_facing };
}
