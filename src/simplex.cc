#include <iostream>
#include <vector>
#include <cfloat>
#include <utils.hh>
#include <vec.hh>
#include <matrix.hh>
using namespace std;

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
                int weWantIdx = 0;
                if (base[i] >= n)
                {
                    for (int j = weWantIdx + 1; j < (int)nonBase.size(); ++j)
                    {
                        if (nonBase[j] < n)
                        {
                            weWantIdx = j;
                            break;
                        }
                    }
                    swap(base[i], nonBase[weWantIdx]);
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
                double tempStep = -solution[base[i]] / db[i];
                if (tempStep < maxStep or (tempStep == maxStep and base[i] < base[exitVariable]))
                {
                    maxStep = tempStep;
                    exitVariable = i;
                }
            }
        }
        /****    MAX STEP: END                ****/
        /**** UPDATES AND BASIS CHANGE: BEGIN ****/
        for (int i = 0; i < (int)base.size(); ++i)
            solution[base[i]] += maxStep * db[i];
        solution[q] = maxStep;
        solution[base[exitVariable]] = 0;

        int idx = 0;
        while (nonBase[idx++] != q)
            continue;
        nonBase[idx] = base[exitVariable];

        base[exitVariable] = q;
        matrix Eta = identityMatrix(m);
        for (int i = 0; i < m; ++i)
        {
            if (i == exitVariable)
                Eta[exitVariable][i] = -1 / db[exitVariable];
            else
                Eta[exitVariable][i] = -db[i] / db[exitVariable];
        }
        Eta = transpose(Eta);
        Binverse = Eta * Binverse;
        z += maxStep * reducedCosts[q];
        /**** UPDATES AND BASIS CHANGE: END   ****/
    }
    return unfeasibleProblem; // decoy
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