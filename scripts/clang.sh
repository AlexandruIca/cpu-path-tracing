#!/usr/bin/env sh

guix environment -m ../scripts/clang-latest.scm --container -- ./clang-build.sh
