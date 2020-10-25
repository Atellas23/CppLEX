#pragma once
#include <vector>
#include <utils.hh>
using namespace std;

using vd = vector<double>;

double operator*(const vd &u, const vd &v);
vd operator+(const vd &u, const vd &v);
vd operator-(const vd &u, const vd &v);
vd operator*(double a, const vd &v);
