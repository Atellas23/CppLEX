#include <matrix.hh>

matrix transpose(const matrix &m)
{
    int n1 = m.size(), n2 = m[0].size();
    matrix res(n2, vd(n1, 0.0));
    for (int i = 0; i < n1; ++i)
    {
        for (int j = 0; j < n2; ++j)
            res[j][i] = m[i][j];
    }
    return res;
}

matrix operator*(const matrix &A, const matrix &B)
{
    int n1 = A.size(), m1 = A[0].size();
    int n2 = B.size(), m2 = B[0].size();
    if (m1 != n2)
    {
        logError(uncoherentArrays);
        return A;
    }
    matrix res(n1, vd(m2, 0));
    for (int i = 0; i < n1; ++i)
    {
        for (int j = 0; j < m2; ++j)
            res[i][j] = A[i] * column(B, j);
    }
    return res;
}

matrix operator*(double a, const matrix &A)
{
    matrix res = A;
    for (int i = 0; i < (int)A.size(); ++i)
    {
        for (int j = 0; j < (int)A[0].size(); ++j)
            res[i][j] *= a;
    }
    return res;
}

matrix takeColumns(const vi &take, const matrix &A)
{
    matrix res;
    for (int i = 0; i < (int)take.size(); ++i)
        res.push_back(column(A, take[i]));
    return transpose(res);
}

matrix identityMatrix(int k)
{
    matrix res(k, vd(k, 0.0));
    for (int i = 0; i < k; ++i)
        res[i][i] = 1;
    return res;
}

vd column(const matrix &A, int idx)
{
    vd res;
    for (vd row : A)
        res.push_back(row[idx]);
    return res;
}

vd operator*(const vd &u, const matrix &B)
{
    int m1 = u.size();
    int n2 = B.size(), m2 = B[0].size();
    if (m1 != n2)
    {
        logError(uncoherentArrays);
        return vd(0);
    }
    vd res(m2, 0.0);

    for (int j = 0; j < m2; ++j)
        res[j] = u * column(B, j);
    return res;
}

vd operator*(const matrix &A, const vd &v)
{
    int n1 = A.size(), m1 = A[0].size();
    int n2 = v.size();
    if (m1 != n2)
    {
        logError(uncoherentArrays);
        return vd(0);
    }
    vd res(n1, 0.0);

    for (int i = 0; i < n1; ++i)
        res[i] = A[i] * v;
    return res;
}

matrix paste(const matrix &A, const matrix &B)
{
    if (A.size() != B.size())
    {
        logError(uncoherentArrays);
        return A;
    }
    matrix res(A.size());
    for (int i = 0; i < (int)A.size(); ++i)
        res[i] = paste(A[i], B[i]);
    return res;
}

void printMatrix(const matrix &A)
{
    for (auto row : A)
        printVec(row);
}