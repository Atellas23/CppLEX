#include <matrix.hh>
#include <iostream>
#include <cmath>
using namespace std;

int main(){
    matrix A = identityMatrix(10);
    matrix expected = {
        {1,0,0,0,0,0,0,0,0,0},
        {0,1,0,0,0,0,0,0,0,0},
        {0,0,1,0,0,0,0,0,0,0},
        {0,0,0,1,0,0,0,0,0,0},
        {0,0,0,0,1,0,0,0,0,0},
        {0,0,0,0,0,1,0,0,0,0},
        {0,0,0,0,0,0,1,0,0,0},
        {0,0,0,0,0,0,0,1,0,0},
        {0,0,0,0,0,0,0,0,1,0},
        {0,0,0,0,0,0,0,0,0,1}
    };

    for (int i = 0; i < (int) expected.size(); i++){
        for (int j = 0; j < (int) expected[0].size(); j++){
            if (abs(A[i][j]-expected[i][j]) > 1e-10){
                cout << "La identidad NO funciona" << endl;
                return 0;
            }
        }
    }
    cout << "La identidad SI funciona" << endl;
}