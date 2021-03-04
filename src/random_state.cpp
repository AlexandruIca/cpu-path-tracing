#include "random_state.hpp"

auto pt::rand_state::default_with_seed(unsigned short const seed) -> rand_state
{
    return rand_state{ std::mt19937{ std::random_device{}() * seed },
                       std::uniform_real_distribution<double>{ 0.0, 1.0 } };
}

auto pt::rand_state::generate() -> double
{
    return dist(rng);
}
