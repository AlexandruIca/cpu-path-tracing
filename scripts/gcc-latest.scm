(use-modules (gnu packages))

(specifications->manifest
  '("cmake"
    "ninja"
    "ccache"
    "gcc-toolchain"))
