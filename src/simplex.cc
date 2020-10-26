#include <iostream>
#include <vector>
#include <cfloat>
#include <utils.hh>
#include <vec.hh>
#include <matrix.hh>
using namespace std;

/* using vd = vector<double>;
using vb = vector<bool>;
using vi = vector<int>;
using matrix = vector<vd>; */

/* enum errorType
{
    unknownError,
    uncoherentArrays
};

enum problemType
{
    unlimitedProblem,
    unfeasibleProblem,
    uniqueOptimum,
    alternativeOptima
};

void logError(errorType err = unknownError)
{
    if (err == uncoherentArrays)
        cerr << "error: non coherent vectors/matrices" << endl;
    else
        cerr << "error: unrecognized operation" << endl;
}*/
/*
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
}*/

/* vd column(const matrix &A, int idx)
{
    vd res;
    for (vd row : A)
        res.push_back(row[idx]);
    return res;
} */

/* matrix transpose(const matrix &m)
{
    int n1 = m.size(), n2 = m[0].size();
    matrix res(n2, vd(n1, 0.0));
    for (int i = 0; i < n1; ++i)
    {
        for (int j = 0; j < n2; ++j)
            res[j][i] = m[i][j];
    }
    return res;
} */

/* vd operator*(const vd &u, const matrix &B)
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
        res[i] = column(A, i) * v;
    return res;
} */

/* matrix operator*(const matrix &A, const matrix &B)
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
} */

/* matrix operator*(int a, const matrix &A)
{
    matrix res(A);
    for (int i = 0; i < (int)A.size(); ++i)
    {
        for (int j = 0; j < (int)A[0].size(); ++j)
            res[i][j] *= a;
    }
    return res;
} */

/* matrix takeColumns(const vi &take, const matrix &A)
{
    matrix res;
    for (int i = 0; i < (int)take.size(); ++i)
        res.push_back(column(A, take[i]));
    return transpose(res);
} */

/* matrix identityMatrix(int k)
{
    matrix res(k, vd(k, 0.0));
    for (int i = 0; i < k; ++i)
        res[i][i] = 1;
    return res;
} */

/* vd paste(const vd &u, const vd &v)
{
    vd res(u);
    for (double d : v)
        res.push_back(d);
    return res;
} */

/* matrix paste(const matrix &A, const matrix &B)
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
} */

/* vd subvec(const vd &u, const vi &idxs)
{
    vd res(idxs.size());
    for (int i = 0; i < (int)idxs.size(); ++i)
        res[i] = u[idxs[i]];
    return res;
} */

/* bool operator>=(const vd &v, double a)
{
    for (double x : v)
    {
        if (x < a)
            return false;
    }
    return true;
} */

/* bool operator>(const vd &v, double a)
{
    for (double x : v)
    {
        if (x <= a)
            return false;
    }
    return true;
} */

int blandRule(const vd &r, const vi &n)
{
    int q = n[0];
    for (int i = 0; i < (int)n.size(); ++i)
    {
        if (r[n[i]] < 0 and n[i] < q)
            q = n[i];
    }
    return q;
}

/* void permutate(vd &x, const vi &perm)
{
    vd backup(x);
    for (int i = 0; i < (int)x.size(); ++i)
        x[i] = backup[perm[i]];
} */

/*
ASP1: PRIMAL SIMPLEX ALGORITHM IMPLEMENTATION
[Using the Bland rule, the product-inverse form, and the initial basic feasible solution (BFS) problem.]

This function takes in a standard form LP problem, and returns its characterization as an unfeasible problem,
an unlimited problem, a unique optimum problem or an alternative optima problem. In any case, it stores in
the reference variable solution the optimum solution it finds, if it exists.

A: matrix of the LP problem in standard form
b: vector of constrictions
n: number of variables in the standard problem
costs: vector of costs of the LP problem
solution: vector where the final solution (if it exists) will be stored
solutionBase: vector where the final solution basis (if it exists) will be stored
*/
problemType ASP1(const matrix &A, const vd &b, int n, const vd &costs, vd &solution, vi &solutionBase, int &iterations, int &phaseIiterations)
{
    int m = A.size(); // rank of A
    /****    INITIALIZATION: BEGIN    ****/

    int phase = 1;

    matrix Ahat = paste(A, identityMatrix(m));

    vd costsHat = paste(vd(n, 0.0), vd(m, 1.0));

    vi base(m);
    for (int i = 0; i < m; ++i)
        base[i] = n + i;

    matrix Binverse = identityMatrix(m);

    vi nonBase(n);
    for (int i = 0; i < n; ++i)
        nonBase[i] = i;

    vd xI = paste(vd(n, 0.0), b); // phase 1 solution candidate

    double z = costsHat * xI;
    /****    INITIALIZATION: END      ****/
    iterations = phaseIiterations = 0;
    cout << "xivato" << endl;
    while (true /*unmet exit conditions*/)
    {
        if (phase == 1)
            ++phaseIiterations;
        ++iterations;
        /****    BFS IDENTIFICATION: BEGIN    ****/
        vd reducedCosts = subvec(costsHat, nonBase) - subvec(costsHat, base) * Binverse * takeColumns(nonBase, A);
        int q = -1;
        if (reducedCosts >= 0)
        {
            if (phase == 2)
                return uniqueOptimum;
            if (z > 0)
                return unfeasibleProblem;

            phase = 2;
            costsHat = costs;

            // check basis inside of {1,...,n}
            for (int i = 0; i < (int)base.size(); ++i)
            {
                if (base[i] >= n)
                {
                    int weWantIdx = 0;
                    for (int j = 0; j < (int)nonBase.size(); ++j)
                    {
                        if (nonBase[j] < n)
                        {
                            weWantIdx = j;
                            break;
                        }
                    }
                    int aux = base[i];
                    base[i] = nonBase[weWantIdx];
                    nonBase[weWantIdx] = aux;
                }
            }
            // and if not, we can pick a non-basic variable and exchange
            // them without taking any effect towards the solution

            // In theory, we have to change the inverse of B as well (?)
            /*
            ||            ||   ||||||||||||||||   ||                      ||
            ||            ||   ||            ||   ||                      ||
            ||            ||   ||            ||   ||                      ||
            ||            ||   ||            ||    \\        //\\        //
            ||||||||||||||||   ||            ||     \\      //  \\      //
            ||            ||   ||            ||      \\    //    \\    //
            ||            ||   ||            ||       \\   ||    ||   //
            ||            ||   ||            ||        \\  //    \\  //
            ||            ||   ||||||||||||||||         \\//      \\//
            */

            Ahat = A;
            z = subvec(costs, base) * subvec(solution, base);
        }
        else
            q = blandRule(reducedCosts, nonBase);
        /****    BFS IDENTIFICATION: END      ****/
        /****    BASIC DIRECTION: BEGIN       ****/
        vd db = (-1) * Binverse * column(A, q);
        if (phase == 2 and db >= 0)
            return unlimitedProblem;
        /****    BASIC DIRECTION: END         ****/
        /****    MAX STEP: BEGIN              ****/
        double maxStep = DBL_MAX;
        int exitVariable = 0;
        for (int i = 0; i < (int)db.size(); ++i)
        {
            if (db[i] < 0)
            {
                if (-solution[base[i]] / db[i] <= maxStep and base[i] < base[exitVariable])
                {
                    maxStep = -solution[base[i]] / db[i];
                    exitVariable = i;
                }
            }
        }
        /****    MAX STEP: END                ****/
        /**** UPDATES AND BASIS CHANGE: BEGIN ****/
        for (int i = 0; i < (int)base.size(); ++i)
            solution[base[i]] += maxStep * db[i];
        solution[q] = maxStep;

        for (int i = 0; i < (int)nonBase.size(); ++i)
        {
            if (nonBase[i] == q)
            {
                nonBase[i] = base[exitVariable];
                break;
            }
        }
        // nonBase[q] = base[exitVariable];
        base[exitVariable] = q;
        matrix Eta = identityMatrix(m);
        vd etarow(m, 0.0);
        for (int i = 0; i < m; ++i)
        {
            if (i == exitVariable)
                etarow[i] = -1 / db[exitVariable];
            etarow[i] = -db[i] / db[exitVariable];
        }
        Eta[exitVariable] = etarow;
        Eta = transpose(Eta);
        Binverse = Eta * Binverse;
        /**** UPDATES AND BASIS CHANGE: END   ****/
    }
    return unfeasibleProblem;
}

int main()
{
    cout << "First implementation of the Simplex algorithm, using" << endl
         << "the Bland rule and the inverse-product form." << endl;
    int n, m;
    cout << "Please enter the dimensions of the matrix (rows, columns): ";
    cin >> n >> m;
    matrix A(n, vd(m, 0.0));
    cout << "Now enter the matrix." << endl;
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < m; ++j)
            cin >> A[i][j];
    }
    cout << "The vector of equality constrictions:" << endl;
    vd b(n, 0.0);
    for (int i = 0; i < n; ++i)
        cin >> b[i];
    cout << "The vector of costs: ";
    vd c(m, 0.0);
    for (int k = 0; k < m; ++k)
        cin >> c[k];
    vd storeSolution(n, 0.0);
    vi storeBase(m, 0);
    int it = 0, phIit = 0;
    cout << "presimplex" << endl;
    problemType PL = ASP1(A, b, n, c, storeSolution, storeBase, it, phIit);
    cout << "postsimplex" << endl;
    if (PL == unlimitedProblem)
    {
        cout << "The problem returned as an unlimited problem after " << phIit << " Simplex phase I iterations and " << it << " phase II iterations." << endl;
        cout << "The found solution is:" << endl;
        for (int i = 0; i < (int)storeSolution.size(); ++i)
            cout << (i ? " " : "") << storeSolution[i];
        cout << endl
             << "for the basis B = {";
        for (int i = 0; i < (int)storeBase.size(); ++i)
            cout << (i ? "," : "") << storeBase[i];
        cout << "}." << endl;
    }
    else if (PL == uniqueOptimum)
    {
        cout << "The problem found an optimum after " << phIit << " Simplex phase I iterations and " << it << " phase II iterations." << endl;
        cout << "The found solution is:" << endl;
        for (int i = 0; i < (int)storeSolution.size(); ++i)
            cout << (i ? " " : "") << storeSolution[i];
        cout << endl
             << "for the basis B = {";
        for (int i = 0; i < (int)storeBase.size(); ++i)
            cout << (i ? "," : "") << storeBase[i];
        cout << "}." << endl;
    }
    else
        cout << "The problem returned as an unfeasible problem after " << phIit << " Simplex phase I iterations." << endl;
}