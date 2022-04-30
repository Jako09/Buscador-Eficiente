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
    file << "f(x)=a*x+b" << '\n';
    file << "g(x)=c*x**2+d*x+e" << '\n';
    file << "FIT_LIMIT=1e-6" << '\n';
    file << "fit f(x) 'Efftime.xxx' u ($1):($3) via a,b" << '\n';
    file << "fit g(x) 'Efftime.xxx' u ($1):($2) via c,d,e" << '\n';
    file << "p 'Efftime.xxx' u ($1):($2) t 'Estandar' lt 1, '' u ($1):($3) t 'Eficiente' lt 2, f(x) t 'Fit Eficiente' lt 2, g(x) t 'Fit Estandar' lt 1" << '\n';
    file << "pause -1" << '\n';
    file.close();
}
