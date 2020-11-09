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
// #define _TEST_READ
// #define CHECK

int read(string filename, matrix &A, vd &b, vd &c)
{
    ifstream file(filename);
    if (not file.good())
        return -1;
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
    return 0;
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
        if (r[i] < 0 and (n[i] < q or rq > 0))
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
problemType ASP1(const matrix &A, const vd &b, const vd &costs, vd &solution, vi &solutionBase, int &iterations, int &phaseIiterations, ofstream &outfile)
{
    outfile.setf(ios::fixed);
    outfile.precision(3);
    outfile << "[CppLEX] Inici de l'ASP amb regla de Bland." << endl;
    outfile << "[CppLEX]\tFase I" << endl;

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
    vd chek, chek2;
    /****    INITIALIZATION: END      ****/
    iterations = phaseIiterations = 0;
    vd reducedCosts, oldSolution = solution;
    while (true /*unmet exit conditions*/)
    {
#ifdef CHECK
        if (oldSolution == solution and iterations != 0)
        {
            cerr << "I've broken at iteration " << ++iterations << endl;
            cerr << "the solution was:" << endl;
            printVec(oldSolution);
            exit(0);
        }
#endif
        if (phase == 1)
            ++phaseIiterations;
        ++iterations;
        /****    BFS IDENTIFICATION: BEGIN    ****/
        reducedCosts = subvec(costsHat, nonBase) - subvec(costsHat, base) * Binverse * takeColumns(nonBase, Ahat);
        int q = -1;
        double rq;
        if (reducedCosts >= 0)
        {
#ifdef CHECK
            // Check if correct
            check = Binverse * takeColumns(base, Ahat);
            assert((int)check.size() == m);
            for (int i = 0; i < m; i++)
                for (int j = 0; j < m; j++)
                    assert((i == j and abs(check[i][j] - 1) < 1e-6) or (i != j and abs(check[i][j]) < 1e-6));
#endif
            if (phase == 2)
            {
                solutionBase = base;
                return uniqueOptimum;
            }
            if (abs(z) > 1e-7)
                return unfeasibleProblem;

            phase = 2;
            costsHat = costs;

            // check basis inside of {1,...,n}
            int weWantIdx = -1;
            for (int i = 0; i < (int)base.size(); ++i)
            {
                if (base[i] >= n) // means base[i] is exiting the base
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
#ifdef CHECK
                    // Check if correct
                    check = Binverse * takeColumns(base, Ahat);
                    assert((int)check.size() == m);
                    for (int i = 0; i < m; i++)
                        for (int j = 0; j < m; j++)
                            assert((i == j and abs(check[i][j] - 1) < 1e-6) or (i != j and abs(check[i][j]) < 1e-6));
#endif
                    swap(base[i], nonBase[weWantIdx]);
                    vd etacolumn(m, 0.0);
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
#ifdef CHECK
                    // Check if correct
                    check = Binverse * takeColumns(base, Ahat);
                    assert((int)check.size() == m);
                    for (int i = 0; i < m; i++)
                        for (int j = 0; j < m; j++)
                            assert((i == j and abs(check[i][j] - 1) < 1e-6) or (i != j and abs(check[i][j]) < 1e-6));
#endif
                }
            }
            while ((int)solution.size() > n)
                solution.pop_back(); // we have to remove artificial variables since we're heading towards phase II

            for (int i = (int)nonBase.size() - 1; i >= 0; --i)
            {
                if (nonBase[i] >= n)
                    nonBase.erase(nonBase.begin() + i);
            }

            // Check solution is valid
            chek = Binverse * bhat;
            chek2 = subvec(solution, base);
            assert(chek.size() == chek2.size());
            for (int i = 0; i < (int)chek.size(); i++)
                assert(abs(chek[i] - chek2[i]) < 1e-5);

            Ahat = A;
            z = subvec(costsHat, base) * subvec(solution, base);
#ifdef CHECK
            // Check if correct
            check = Binverse * takeColumns(base, Ahat);
            assert((int)check.size() == m);
            for (int i = 0; i < m; i++)
                for (int j = 0; j < m; j++)
                    assert((i == j and abs(check[i][j] - 1) < 1e-6) or (i != j and abs(check[i][j]) < 1e-6));
#endif
            outfile << "[CppLEX]\t\tSBF inicial trobada a la iteració " << iterations << endl
                    << "[CppLEX]\tFase II" << endl;
            continue;
        }
        else
        {
            q = blandRule(reducedCosts, nonBase, rq);
#ifdef CHECK
            assert(rq < 0);
#endif
        }
#ifdef CHECK
        // Check solution is valid
        chek = Binverse * bhat;
        chek2 = subvec(solution, base);
        assert(chek.size() == chek2.size());
        for (int i = 0; i < (int)chek.size(); i++)
            assert(abs(chek[i] - chek2[i]) < 1e-5);
#endif
        /****    BFS IDENTIFICATION: END      ****/
        /****    BASIC DIRECTION: BEGIN       ****/
        vd db = (-1) * Binverse * column(Ahat, q);
        if (phase == 2 and db >= 0)
        {
            solutionBase = base;
            return unlimitedProblem;
        }
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

        outfile << "[CppLEX]\t\titer " << iterations << ": q = " << q << ", rq = " << rq << ", B(p) = " << base[exitVariable] << ", \u03B8* = " << maxStep;

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
#ifdef CHECK
        // Check if correct
        check = Binverse * takeColumns(base, Ahat);
        assert((int)check.size() == m);
        for (int i = 0; i < m; i++)
            for (int j = 0; j < m; j++)
                assert((i == j and abs(check[i][j] - 1) < 1e-6) or (i != j and abs(check[i][j]) < 1e-6));
#endif
        z += maxStep * rq;

        outfile << ", z=" << z << endl;
#ifdef CHECK
        // Check solution is valid
        chek = Binverse * bhat;
        chek2 = subvec(solution, base);
        assert(chek.size() == chek2.size());
        for (int i = 0; i < (int)chek.size(); i++)
            assert(abs(chek[i] - chek2[i]) < 1e-5);
#endif
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
    matrix A;
    vd b, c;
    int err = read(argv[1], A, b, c);
    if (err < 0)
    {
        cerr << "error: file \"" << argv[1] << "\" could not be loaded" << endl;
        exit(0);
    }
#ifndef _TEST_READ
    vd storeSolution(c.size(), 0.0);
    vi storeBase(b.size(), 0);
    int it = 0, phIit = 0;
    ofstream outputDataFile(string(argv[1]) + ".out");
    problemType PL = ASP1(A, b, c, storeSolution, storeBase, it, phIit, outputDataFile);
    if (PL == unlimitedProblem)
    {
        outputDataFile << "[CppLEX]\t\tHem determinat que el problema és il·limitat a la iteració " << it << endl;
        outputDataFile << "[CppLEX] Fi de l'ASP" << endl
                       << endl
                       << endl
                       << "B* = " << printVec(storeBase) << endl
                       << "xB* = " << printVec(subvec(storeSolution, storeBase)) << endl;
        cout << "The problem returned as an unlimited problem after " << phIit << " Simplex phase I iterations and " << it << " phase I + phase II iterations." << endl;
        cout << "The found solution is:" << endl
             << printVec(subvec(storeSolution, storeBase)) << endl;
        cout << endl
             << "for the basis B* = " << printVec(storeBase) << endl;
    }
    else if (PL == uniqueOptimum)
    {
        outputDataFile << "[CppLEX]\t\tSolució òptima trobada, iteració " << it << ", z=" << c * storeSolution << endl;
        outputDataFile << "[CppLEX] Fi de l'ASP" << endl
                       << endl
                       << endl
                       << "B* = " << printVec(storeBase) << endl
                       << "xB* = " << printVec(subvec(storeSolution, storeBase)) << endl;
        cout << "The problem found an optimum after " << phIit << " Simplex phase I iterations and " << it << " phase II iterations." << endl;
        cout << "The found solution is xB* = " << endl
             << printVec(subvec(storeSolution, storeBase));
        cout << endl
             << "for the basis B* = " << printVec(storeBase) << endl;
    }
    else
    {
        cout << "The problem returned as an unfeasible problem after " << phIit << " Simplex phase I iterations." << endl;
        outputDataFile << "[CppLEX]\t\tHem determinat que el problema és infactible a la iteració " << it << endl;
        outputDataFile << "[CppLEX] Fi de l'ASP" << endl;
    }
#endif
}