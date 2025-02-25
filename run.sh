#!/bin/bash

g++ -std=c++17 main.cpp sources/amp.cpp sources/utils.cpp -Ofast -o main -lcuba -lm