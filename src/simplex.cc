#include <iostream>
#include <vector>
#include <cfloat>
#include <fstream>
#include <utils.hh>
#include <vec.hh>
#include <matrix.hh>
#include <cassert>
#include <cmath>
using namespace std;
//#define _TEST_READ

void read(string filename, matrix &A, vd &b, vd &c)
{
    ifstream file(filename);
    if (not file.good())
    {
        cerr << "error: the file \"" << filename << "\" could not be loaded." << endl;
        return;
    }
    bool matr = true, restrictions = false, costs = false;
    vd tempVec;
    while (not file.eof())
    {
        string buff;
        getline(file, buff);
        trimWhitespace(buff);
        if (buff == "end of matrix")
        {
            matr = false;
            restrictions = true;
            continue;
        }
        else if (buff == "end of restrictions")
        {
            restrictions = false;
            costs = true;
            continue;
        }
        else if (buff == "end of costs")
        {
            costs = false;
            break;
        }
        vector<string> tkn = tokenize(buff);
        if (matr)
        {
            for (auto a : tkn)
                tempVec.push_back(stod(a));
            A.push_back(tempVec);
            tempVec.clear();
        }
        else if (restrictions)
            for (auto bi : tkn)
                b.push_back(stod(bi));
        else if (costs)
            for (auto ci : tkn)
                c.push_back(stod(ci));
    }
#ifdef _TEST_READ
    cout << "Matrix A:" << endl;
    for (auto row : A)
    {
        for (int i = 0; i < (int)row.size(); ++i)
            cout << (i ? " " : "") << row[i];
        cout << endl;
    }
    cout << endl
         << "Restrictions b:" << endl;
    for (int j = 0; j < (int)b.size(); ++j)
        cout << (j ? " " : "") << b[j];
    cout << endl
         << endl
         << "Costs c:" << endl;
    for (int j = 0; j < (int)c.size(); ++j)
        cout << (j ? " " : "") << c[j];
    cout << endl;
#endif
}

// Precondition: There is one negative value in r
int blandRule(const vd &r, const vi &n, double &rq)
{
    int q = n[0];
    rq = r[0];
    for (int i = 1; i < (int)n.size(); ++i)
    {
        if (r[i] < 0 and n[i] < q)
        {
            q = n[i];
            rq = r[i];
        }
    }
    return q;
}

// Precondition: There is one negative value in r
int minValueRule(const vd &r, const vi &n)
{
    int q = n[0];
    for (int i = 1; i < (int)n.size(); i++)
    {
        if (r[n[i]] < r[q])
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
b: non-negative vector of constrictions
n: number of variables in the standard problem
costs: vector of costs of the LP problem
solution: vector where the final solution (if it exists) will be stored
solutionBase: vector where the final solution basis (if it exists) will be stored
*/
problemType ASP1(const matrix &A, const vd &b, const vd &costs, vd &solution, vi &solutionBase, int &iterations, int &phaseIiterations)
{
    int n = costs.size();
    int m = b.size(); // rank of A
    /****    INITIALIZATION: BEGIN    ****/

    int phase = 1;

    matrix Ahat = paste(A, identityMatrix(m));

    vd costsHat = paste(vd(n, 0.0), vd(m, 1.0));

    vd bhat = b;
    // checking that b is non-negative (else phase I does not finish)
    for (int i = 0; i < m; ++i)
    {
        if (b[i] < 0)
        {
            bhat[i] *= -1;
            Ahat[i] = (-1) * Ahat[i]; // Ax is EQUAL to b so we don't have to care for changing the polarity of inequalities
        }
    }

    vi base(m);
    for (int i = 0; i < m; ++i)
        base[i] = n + i;

    matrix Binverse = identityMatrix(m);

    vi nonBase(n);
    for (int i = 0; i < n; ++i)
        nonBase[i] = i;

    solution = paste(vd(n, 0.0), bhat); // phase 1 solution candidate

    double z = costsHat * solution;
    matrix check;
    /****    INITIALIZATION: END      ****/
    iterations = phaseIiterations = 0;
    vd reducedCosts, oldSolution = solution;
    while (true /*unmet exit conditions*/)
    {
        if (oldSolution == solution and iterations != 0)
        {
            cerr << "I've broken at iteration " << ++iterations << endl;
            cerr << "the solution was:" << endl;
            printVec(oldSolution);
            exit(0);
        }
        if (phase == 1)
            ++phaseIiterations;
        ++iterations;
        /****    BFS IDENTIFICATION: BEGIN    ****/
        // cout << "[iter " << iterations << "] Matrix A:" << endl;
        //printMatrix(Ahat);
        // cout << "[iter " << iterations << "] Matrix B^(-1):" << endl;
        //printMatrix(Binverse);
        cout << "[iter " << iterations << "] Base:" << endl;
        printVec(base);
        // cout << "[iter " << iterations << "] nonBase:" << endl;
        // printVec(nonBase);
        reducedCosts = subvec(costsHat, nonBase) - subvec(costsHat, base) * Binverse * takeColumns(nonBase, Ahat);
        cout << "[iter " << iterations << "] reducedCosts:" << endl;
        printVec(reducedCosts);
        // cout << "pot calcular reducedCosts" << endl;
        cout << "z: " << z << endl;
        int q = -1;
        double rq;
        if (reducedCosts >= 0)
        {
            // Check if correct
            check = Binverse * takeColumns(base, Ahat);
            assert((int)check.size() == m);
            for (int i = 0; i < m; i++)
                for (int j = 0; j < m; j++)
                    assert((i == j and abs(check[i][j] - 1) < 1e-6) or (i != j and abs(check[i][j]) < 1e-6));

            if (phase == 2)
                return uniqueOptimum;
            if (z > 0)
                return unfeasibleProblem;

            phase = 2;
            costsHat = costs;
            // cout << "hem entrat a fase II" << endl;

            // check basis inside of {1,...,n}
            int weWantIdx = -1;
            for (int i = 0; i < (int)base.size(); ++i)
            {
                if (base[i] >= m) // means base[i] is exiting the base
                {
                    for (int j = weWantIdx + 1; j < (int)nonBase.size(); ++j)
                    {
                        if (nonBase[j] < n)
                        {
                            weWantIdx = j; // the variable in nonBase[weWantIdx] is the one that will enter the basis
                            break;
                        }
                    }
                    /* CHANGING THE BASIC MATRIX WHEN EXCHANGING DEGENERATE ORIGINAL VARIABLES: START */
                    vd dullBasicDirection = (-1) * Binverse * column(Ahat, nonBase[weWantIdx]);
                    // here, maxStep is 0 because we're not really moving between extreme points (we have degeneration)
                    solution[nonBase[weWantIdx]] = solution[base[i]] = 0;

                    // Check if correct
                    check = Binverse * takeColumns(base, Ahat);
                    assert((int)check.size() == m);
                    for (int i = 0; i < m; i++)
                        for (int j = 0; j < m; j++)
                            assert((i == j and abs(check[i][j] - 1) < 1e-6) or (i != j and abs(check[i][j]) < 1e-6));

                    swap(base[i], nonBase[weWantIdx]);
                    vd etacolumn(m, 0.0);
                    // exitVariable = i
                    for (int k = 0; k < m; ++k)
                    {
                        if (k == i)
                            etacolumn[k] = -1 / dullBasicDirection[i];
                        else
                            etacolumn[k] = -dullBasicDirection[k] / dullBasicDirection[i];
                    }
                    vd aux = Binverse[i];
                    for (int l = 0; l < m; ++l)
                    {
                        if (l != i)
                            Binverse[l] = Binverse[l] + etacolumn[l] * aux;
                        else
                            Binverse[l] = etacolumn[l] * Binverse[l];
                    }
                    /* CHANGING THE BASIC MATRIX WHEN EXCHANGING DEGENERATE ORIGINAL VARIABLES: END */
                    // Check if correct
                    check = Binverse * takeColumns(base, Ahat);
                    assert((int)check.size() == m);
                    for (int i = 0; i < m; i++)
                        for (int j = 0; j < m; j++)
                            assert((i == j and abs(check[i][j] - 1) < 1e-6) or (i != j and abs(check[i][j]) < 1e-6));
                }
            }
            while ((int)solution.size() > n)
                solution.pop_back(); // we have to remove artificial variables since we're heading towards phase II

            for (int i = (int)nonBase.size() - 1; i >= 0; --i)
            {
                if (nonBase[i] >= n)
                    nonBase.erase(nonBase.begin() + i);
            }

            Ahat = A;
            z = subvec(costs, base) * subvec(solution, base);
            // Check if correct
            check = Binverse * takeColumns(base, Ahat);
            assert((int)check.size() == m);
            for (int i = 0; i < m; i++)
                for (int j = 0; j < m; j++)
                    assert((i == j and abs(check[i][j] - 1) < 1e-6) or (i != j and abs(check[i][j]) < 1e-6));
            cout << "Entered Phase II correctly" << endl;
            continue;
        }
        else
        {
            // cout << "apliquem la regla de bland" << endl;
            q = blandRule(reducedCosts, nonBase, rq);
            assert(rq < 0);
            // cout << "rq: " << rq << endl;
            // cout << "hem acabat d'aplicar la regla de bland" << endl;
        }
        /****    BFS IDENTIFICATION: END      ****/
        /****    BASIC DIRECTION: BEGIN       ****/
        // cout << "[iter " << iterations << "] Matrix Ahat:" << endl;
        //printMatrix(Ahat);
        vd db = (-1) * Binverse * column(Ahat, q);
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
                /* cout << "BASE:" << endl;
                printVec(base);
                cout << endl
                     << "DIRECCIO BASICA:" << endl;
                printVec(db);
                cout << endl
                     << "SOLUCIO:" << endl;
                printVec(solution); */
                double tempStep = -solution[base[i]] / db[i];
                if (tempStep < maxStep or (tempStep == maxStep and base[i] < base[exitVariable]))
                {
                    maxStep = tempStep;
                    exitVariable = i;
                }
            }
        }
        cout << "\u03B8*=" << maxStep << endl;
        cout << "B(p)=" << base[exitVariable] << endl;
        cout << "q=" << q << endl;
        cout << "[iter " << iterations << "] db:" << endl;
        printVec(db);
        cout << "CURRENT SOLUTION (" << solution.size() << "): ";
        printVec(solution);
        /****    MAX STEP: END                ****/
        /**** UPDATES AND BASIS CHANGE: BEGIN ****/
        oldSolution = solution;
        for (int i = 0; i < (int)base.size(); ++i)
            solution[base[i]] += maxStep * db[i];
        solution[q] = maxStep;
        solution[base[exitVariable]] = 0;

        int idx = 0;
        while (nonBase[idx++] != q)
            continue;
        nonBase[--idx] = base[exitVariable];

        base[exitVariable] = q;

        vd etacolumn(m, 0.0);
        cout << "exitVariable: " << exitVariable << endl;
        for (int i = 0; i < m; ++i)
        {
            if (i == exitVariable)
                etacolumn[i] = -1 / db[exitVariable];
            else
                etacolumn[i] = -db[i] / db[exitVariable];
        }
        vd aux = Binverse[exitVariable];
        for (int i = 0; i < m; ++i)
        {
            if (i != exitVariable)
                Binverse[i] = Binverse[i] + etacolumn[i] * aux;
            else
                Binverse[i] = etacolumn[i] * Binverse[i];
        }
        // Check if correct
        check = Binverse * takeColumns(base, Ahat);
        assert((int)check.size() == m);
        for (int i = 0; i < m; i++)
            for (int j = 0; j < m; j++)
                assert((i == j and abs(check[i][j] - 1) < 1e-6) or (i != j and abs(check[i][j]) < 1e-6));

        z += maxStep * rq;
        /**** UPDATES AND BASIS CHANGE: END   ****/
    }
    return unfeasibleProblem; // decoy
}

int main(int argc, char **argv)
{
    if (argc < 2)
    {
        cerr << "usage: " << argv[0] << " <data file name>" << endl;
        exit(0);
    }
    /* cout << "First implementation of the Simplex algorithm, using" << endl
         << "the Bland rule and the inverse-product form." << endl; */
    matrix A;
    vd b, c;
    read(argv[1], A, b, c);
#ifndef _TEST_READ
    vd storeSolution(c.size(), 0.0);
    vi storeBase(b.size(), 0);
    int it = 0, phIit = 0;
    cout << "presimplex" << endl;
    problemType PL = ASP1(A, b, c, storeSolution, storeBase, it, phIit);
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
#endif
}