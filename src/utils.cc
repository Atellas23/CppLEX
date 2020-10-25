#include <utils.hh>

void logError(errorType err = unknownError)
{
    if (err == uncoherentArrays)
        std::cerr << "error: non coherent vectors/matrices" << std::endl;
    else
        std::cerr << "error: unrecognized operation" << std::endl;
}