#!/bin/bash
g++ -o keisan main.cpp hamiltonian.cpp kubostreda.cpp matrix.cpp -fopenmp -llapack -lblas -lm
