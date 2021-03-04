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

    [[nodiscard]] static auto default_with_seed(unsigned short const seed) -> rand_state
    {
        return rand_state{ std::mt19937{ std::random_device{}() * seed },
                           std::uniform_real_distribution<double>{ 0.0, 1.0 } };
    }

    [[nodiscard]] auto generate() -> double
    {
        return dist(rng);
    }
};

} // namespace pt

#endif // !PT_RAND_STATE
