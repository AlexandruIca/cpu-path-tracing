#!/usr/bin/env sh

# Invoke this script with guix:
# guix environment -m ../scripts/gcc-latest.scm --container -- ./gcc-build.sh

CC=clang CXX=clang++ cmake -G"Ninja" ..

mv compile_commands.json ..

cmake --build .
