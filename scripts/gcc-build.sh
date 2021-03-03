#!/usr/bin/env sh

# Invoke this script with guix:
# guix environment -m ../scripts/gcc-latest.scm --container -- ./gcc-build.sh

CC=gcc CXX=g++ cmake -G"Ninja" ..

cmake --build .
