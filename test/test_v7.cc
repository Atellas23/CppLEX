#include <vec.hh>
#include <iostream>
#include <cmath>
using namespace std;

int main(){
    vd u = {0.106652770180584,0.961898080855054,0.004634224134067,0.774910464711502,0.817303220653433,0.868694705363510,0.084435845510910,0.399782649098896,0.259870402850654,0.800068480224308};
    double a = 0, b = 1, c = 0.5, d = 0.004634224134067;
    if (u >= b or u >= c or not (u >= a) or not (u >= d))
        cout << "El >= NO funciona" << endl;
    else
        cout << "El >= SI funciona" << endl;
}