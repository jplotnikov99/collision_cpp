#!/bin/bash

g++ -std=c++17 main.cpp sources/amp.cpp sources/ampmassless.cpp sources/utils.cpp -Ofast -march=native -ffast-math -funroll-loops -fomit-frame-pointer -DNDEBUG -fopenmp -o main -lcuba -lm