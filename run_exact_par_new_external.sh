#!/bin/bash
g++ -O3 -std=c++11 -lgomp -fopenmp /home/isis/MoCHy/main_exact_par2.cpp -o /home/isis/MoCHy/exact_par2;
/home/isis/MoCHy/exact_par2 64 $1 $2;
rm /home/isis/MoCHy/exact_par2;
