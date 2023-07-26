# g++-12 -O3 -std=c++11 -lgomp -fopenmp main_exact_par2.cpp -o exact_par;
g++ -O3 -std=c++11 -lgomp -fopenmp main_exact_par2.cpp -o exact_par;
# ./exact_par 8 $1 $2;
./exact_par 4;
rm exact_par;
