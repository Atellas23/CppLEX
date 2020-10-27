#include <vec.hh>
#include <iostream>
#include <cmath>
using namespace std;

int main(){
    vd u = {0.106652770180584,0.961898080855054,0.004634224134067,0.774910464711502,0.817303220653433,0.868694705363510,0.084435845510910,0.399782649098896,0.259870402850654,0.800068480224308};
    double a = 0.853031117721894;
    vd expected = {0.090978131755280,0.820528995046331,0.003953137392857,0.661022739847245, 0.697185079831701,0.741023615475326,0.072026403671965,0.341027040006651,0.221677540206532,0.682483309939798};
    
    vd prod = a*u;
    for (int i = 0; i < (int) u.size(); i++)
        if (abs(prod[i] - expected[i]) > 1e-10){
            cout << "El producto por escalar NO funciona" << endl;
            return 0;     
        }
    cout << "El producto por escalar SI funciona" << endl;
}