#pragma once
#include <iostream>
#include <string>
#include <sstream>
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
    uniqueOptimum
};

void logError(errorType err);
void trimWhitespace(string &s);
vector<string> tokenize(string &s);