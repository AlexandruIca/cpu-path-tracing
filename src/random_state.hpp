#ifndef PT_RAND_STATE
#define PT_RAND_STATE
#pragma once

#include <algorithm>
#include <array>
#include <functional>
#include <random>

namespace pt {

struct rand_state
{
    std::mt19937 rng;
    std::uniform_real_distribution<double> dist;

    rand_state() = delete;

    [[nodiscard]] static auto default_with_seed(unsigned short seed) -> rand_state;
    [[nodiscard]] auto generate() -> double;
};

} // namespace pt

#endif // !PT_RAND_STATE
