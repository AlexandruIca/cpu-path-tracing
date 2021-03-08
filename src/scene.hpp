#ifndef PT_SCENE_HPP
#define PT_SCENE_HPP
#pragma once

#include <vector>

#include "camera.hpp"
#include "sphere.hpp"

namespace pt {

struct scene
{
    std::vector<sphere> spheres{};
    camera_config camera_parameters{};
};

} // namespace pt

#endif // !PT_SCENE_HPP
