#!/bin/bash
g++ -O3 -std=c++11 -lgomp -fopenmp /home/isis/MoCHy/main_exact_par.cpp -o /home/isis/MoCHy/exact_par;
/home/isis/MoCHy/exact_par 32;
rm /home/isis/MoCHy/exact_par;