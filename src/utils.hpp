#ifndef PT_UTILS_HPP
#define PT_UTILS_HPP
#pragma once

namespace pt {

///
/// \brief Clamps x to (0.0, 1.0).
///
[[nodiscard]] auto clamp(double x) noexcept -> double;
///
/// \brief Returns gamma 2.2 corrected color in range 0-256.
///
[[nodiscard]] auto color_to_int(double x) noexcept -> int;

} // namespace pt

#endif // !PT_UTILS_HPP
