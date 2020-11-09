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

vd paste(const vd &u, const vd &v)
{
    vd res(u);
    for (double d : v)
        res.push_back(d);
    return res;
}

vd subvec(const vd &u, const vi &idxs)
{
    vd res(idxs.size());
    for (int i = 0; i < (int)idxs.size(); ++i)
        res[i] = u[idxs[i]];
    return res;
}

bool operator>=(const vd &v, double a)
{
    for (double x : v)
    {
        if (x < a)
            return false;
    }
    return true;
}

bool operator>(const vd &v, double a)
{
    for (double x : v)
    {
        if (x <= a)
            return false;
    }
    return true;
}

string printVec(const vd &a)
{
    stringstream ss;
    ss << "[";
    for (int i = 0; i < (int)a.size(); ++i)
        ss << (i ? ", " : "") << a[i];
    ss << "]";
    return ss.str();
}

string printVec(const vi &a)
{
    stringstream ss;
    ss << "{";
    for (int i = 0; i < (int)a.size(); ++i)
        ss << (i ? ", " : "") << a[i];
    ss << "}";
    return ss.str();
}