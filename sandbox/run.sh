#!/usr/bin/env sh

g++ -O3 -fopenmp main.cpp -o smallpt

time ./smallpt $1
