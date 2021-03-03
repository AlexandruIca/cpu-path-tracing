#!/usr/bin/env sh

CLANG_FORMAT_PATH=clang-format ./scripts/clang-format.sh

./scripts/check-cmake-format.sh

CLANG_TIDY_PATH=clang-tidy ./scripts/clang-tidy.sh
