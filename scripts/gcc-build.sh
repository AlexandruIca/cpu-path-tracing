#!/usr/bin/env sh

CC=gcc CXX=g++ cmake -G"Ninja" \
    -DENABLE_SANITIZER_ADDRESS=ON \
    -DENABLE_SANITIZER_UNDEFINED_BEHAVIOR=ON \
    ..

cmake --build .
