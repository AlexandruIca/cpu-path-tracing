(use-modules (gnu packages))

(specifications->manifest
  '("cmake"
    "ninja"
    "ccache"
    "libomp"
    "clang-toolchain"))
