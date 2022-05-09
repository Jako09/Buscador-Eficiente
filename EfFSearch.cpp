// Implementation for efficient neighbor finder
// Author: Jan Carlo Alvarez Centeno
#include"Dis0.h"
#include<chrono>
#include"gnuplot.h"
using namespace std;

int main(){
    //There is four steps to do the efficient neighbor finder
    // 1.- Discretice the space in cells and give each cell a unique key
    // 2.- Order the keys obteined by each particle from least to greatest
    // 3.- Classify each particle according to key
    // 4.- Search in cells according to the smoothing lenght h
    int Ptype=2; // Type Problem, 1.- is for partilce in a box; 2.- is for Harmonic Oscillator
    int N, N_itn, nclass, xf;
    int t, tEff; // for time
    string formatimage, term;
    printf("¿Cuántas iteraciones deseas realizar? Comenzamos desde 100 y las iteraciones van en multiplos de 100 \n");
    cin >> N_itn;
    printf("Introduce el valor inicial del número de partículas N \n");
    cin >> N;
    printf("¿Qué tipo de archivo de salida deseas? \n EPS\t ó \tPNG, escribe la terminación en mayusculas \n");
    cin >> formatimage;
    printf("¿Qué sistema operativo usas?\n Si es windows, escribe 'windows'\n Si es linux, escribe 'x11'\n");
    cin >> term;
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::microseconds;
    double h, h_hash, xmin, xmax;
    xmax=4.0;
    xmin=-4.0;
    xf=1;
    //here will be the function to calculate the xf, where xf is number of neight cells that will check for each particle
    ofstream file("Efftime.dat");
    ofstream file1("Densidad.dat");
    // First we need to know the position of the variables, so we call any function to initialize the distribution of the particles
    for(int i=1; i<N_itn;i++){
        N+=100;
        double m[N]={}, R[N]={}, D[N]={}, DEff[N]={};
        int keyS[N]={}, idx[N]={};
        grid(Ptype, N, m, R);
        h=200.0/(double)N;
        h_hash=2*h; //Only for the standard searching
        nclass=int((xmax-xmin)/h_hash);
        int idxmin[nclass]={}, idxmax[nclass]={};
        int act[nclass]={};
        auto start = high_resolution_clock::now();
        Densidad0(N, m, R, h, D); //In general the order is the dependency of the variables, by the adapptative smoothing lenght exist a dependency for R
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        t=duration.count(); //time in microseconds
        auto startEff = high_resolution_clock::now();
        CensoSPH(N, h_hash, xmin, R, keyS, idx, idxmin,  idxmax, act);
        Densidad0Eff( N,nclass, xf, m, R, h,  DEff, keyS,idx, idxmin, idxmax, act);
        auto stopEff = high_resolution_clock::now();
        auto durationEff = duration_cast<microseconds>(stopEff - startEff);
        tEff=durationEff.count(); //time in microseconds
        //-------------This is for time vs #particles print
        file << N << " " << t << " " << tEff << '\n';

        //-------------This is for density print
        file1 << "\n\n\n";
        for(int i=0; i<N; i++){
            file1 << R[i] << "\t" << D[i] << "\t" <<  DEff[i] << "\n";
        }
        printf("iteracion %d \n ", i);

    }
    printGraphDensity(formatimage, N_itn, term);
    printGraphtimeEff(formatimage);
    file.close();
    file1.close();
    cout << "The process finished" << '\n';
    return 0;
}
