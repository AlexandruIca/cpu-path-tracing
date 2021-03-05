#include "sphere.hpp"
#include "constants.hpp"

#include <cmath>

auto pt::sphere::intersect(ray const& r) const noexcept -> double
{
    vec3 op = position - r.origin;
    double t = 0.0;
    double const b = op.dot(r.direction);
    double det = b * b - op.dot(op) + radius * radius;

    if(det < 0) {
        return 0;
    }

    det = std::sqrt(det);
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
