#!/usr/bin/env sh

CC=gcc CXX=g++ cmake -G"Ninja" \
    -DCMAKE_BUILD_TYPE=Release \
    ..

cmake --build .
