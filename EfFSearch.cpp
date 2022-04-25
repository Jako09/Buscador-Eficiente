// Implementation for efficient neighbor finder
// Author: Jan Carlo Alvarez Centeno
#include"Dis0.h"
#include<chrono>
using namespace std;

int main(){
    //There is four steps to do the efficient neighbor finder
    // 1.- Discretice the space in cells and give each cell a unique key
    // 2.- Order the keys obteined by each particle from least to greatest
    // 3.- Classify each particle according to key
    // 4.- Search in cells according to the smoothing lenght h
    int Ptype=2; // Type Problem, 1.- is for partilce in a box; 2.- is for Harmonic Oscillator
    int N=640;
    int t; // for time
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::microseconds;
    double m[N], h[N], R[N], D[N];
    // First we need to know the position of the variables, so we call any function to initialize the distribution of the particles
    grid(Ptype, N, m, R);
    for(int i=0; i<N; i++){
        h[i]=200.0/(double)N;
    }
    auto start = high_resolution_clock::now();
    Densidad0(N, m, R, h, D); //In general the order is the dependency of the variables, by the adapptative smoothing lenght exist a dependency for R
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    t=duration.count(); //time in microseconds
    printf("The time duration for the Densidad function it was: %d μs \n",t);
}
