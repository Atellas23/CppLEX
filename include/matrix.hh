#pragma once
#include <vec.hh>
#include <utils.hh>
using namespace std;

using matrix = vector<vd>;

matrix transpose(const matrix &m);
matrix operator*(const matrix &A, const matrix &B);
matrix operator*(double a, const matrix &A);
matrix takeColumns(const vi &take, const matrix &A);
matrix identityMatrix(int k);
vd column(const matrix &A, int idx);
vd operator*(const vd &u, const matrix &B);
vd operator*(const matrix &A, const vd &v);
matrix paste(const matrix &A, const matrix &B);
void printMatrix(const matrix &A);