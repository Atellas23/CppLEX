#include <vec.hh>
#include <iostream>
#include <cmath>
using namespace std;

int main(){
    vd u = {0.106652770180584,0.961898080855054,0.004634224134067,0.774910464711502,0.817303220653433,0.868694705363510,0.084435845510910,0.399782649098896,0.259870402850654,0.800068480224308};
    vi idx = {0,4,2,3};
    vd expected = {0.106652770180584,0.817303220653433,0.004634224134067,0.774910464711502};
    
    vd index = subvec(u, idx);
    for (int i = 0; i < (int) index.size(); i++)
        if (abs(index[i] - expected[i]) > 1e-10){
            cout << "Subvec NO funciona" << endl;
            return 0;     
        }
    cout << "Subvec SI funciona" << endl;
}