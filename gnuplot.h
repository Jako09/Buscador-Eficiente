#include<fstream>

using namespace std;

void printGraphEff(){
    ofstream file("SearchEff.gnu");
    file << "set term postscript eps enhanced color" << '\n';
    file << "set out 'SearchEff.eps'" << '\n';
    file << "set grid" << '\n';
    file << "set title 'Tiempo vs Numero de particulas'" << '\n';
    file << "set ylabel 't:= tiempo de computo'"<<'\n';
    file << "set xlabel 'N:= numero de particulas'"<<'\n';
    file << "p 'Efftime.xxx' u ($1):($3), '' u ($2):($3)" << '\n';
    file << "pause -1" << '\n';
    file.close();
}