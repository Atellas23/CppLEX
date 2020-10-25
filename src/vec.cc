#include <vec.hh>

double operator*(const vd &u, const vd &v)
{
    if (u.size() != v.size())
    {
        logError(uncoherentArrays);
        return 0;
    }
    double s = 0;
    for (int i = 0; i < (int)u.size(); ++i)
        s += u[i] * v[i];
    return s;
}

vd operator+(const vd &u, const vd &v)
{
    if (u.size() != v.size())
    {
        logError(uncoherentArrays);
        return u;
    }
    vd res(u);
    for (int i = 0; i < (int)u.size(); ++i)
        res[i] += v[i];
    return res;
}

vd operator-(const vd &u, const vd &v)
{
    if (u.size() != v.size())
    {
        logError(uncoherentArrays);
        return u;
    }
    vd res(u);
    for (int i = 0; i < (int)u.size(); ++i)
        res[i] -= v[i];
    return res;
}

vd operator*(double a, const vd &v)
{
    vd res(v);
    for (double &x : res)
        x *= a;
    return res;
}