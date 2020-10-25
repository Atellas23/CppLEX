#pragma once
#include <iostream>

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
