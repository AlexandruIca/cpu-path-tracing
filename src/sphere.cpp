#include "sphere.hpp"
#include "constants.hpp"

#include <cmath>

auto pt::sphere::intersect(ray const& r) const noexcept -> double
{
    vec3 const oc = r.origin - position;
    auto const a = r.direction.dot(r.direction);
    auto const half_b = oc.dot(r.direction);
    auto const c = oc.dot(oc) - radius * radius;
    auto const discriminant = half_b * half_b - a * c;

    if(discriminant < 0) {
        return 0.0;
    }

    auto const sqrtd = std::sqrt(discriminant);
    auto root = (-half_b - sqrtd) / a;

    if(root < epsilon) {
        root = (-half_b + sqrtd) / a;

        if(root < epsilon) {
            return 0.0;
        }
    }

    return root;
}
