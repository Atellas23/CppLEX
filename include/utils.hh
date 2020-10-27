#pragma once
#include <iostream>
#include <string>
#include <vector>
using namespace std;

enum errorType
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

void logError(errorType err);
void trimWhitespace(string &s);
vector<string> tokenize(string &s);