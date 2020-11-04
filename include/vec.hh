#pragma once
#include <vector>
#include <utils.hh>
using namespace std;

using vd = vector<double>;
using vb = vector<bool>;
using vi = vector<int>;

double operator*(const vd &u, const vd &v);
vd operator+(const vd &u, const vd &v);
vd operator-(const vd &u, const vd &v);
vd operator*(double a, const vd &v);
vd paste(const vd &u, const vd &v);
vd subvec(const vd &u, const vi &idxs);
bool operator>=(const vd &v, double a);
bool operator>(const vd &v, double a);
void printVec(const vd &a);
void printVec(const vi &a);
