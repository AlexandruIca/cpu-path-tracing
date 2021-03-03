#!/usr/bin/env sh

CC=clang CXX=clang++ cmake -G"Ninja" \
    -DENABLE_SANITIZER_ADDRESS=ON \
    -DENABLE_SANITIZER_UNDEFINED_BEHAVIOR=ON \
    ..

mv compile_commands.json ..

cmake --build .
