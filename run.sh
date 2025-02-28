#!/bin/bash

g++ -std=c++17 main.cpp sources/amp.cpp sources/ampmassless.cpp sources/utils.cpp -Ofast -o main -lcuba -lm