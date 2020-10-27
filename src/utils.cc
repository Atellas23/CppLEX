#include <utils.hh>
using namespace std;

void logError(errorType err = unknownError)
{
    if (err == uncoherentArrays)
        std::cerr << "error: non coherent vectors/matrices" << std::endl;
    else
        std::cerr << "error: unrecognized operation" << std::endl;
}

void trimWhitespace(string &s)
{
    while (s[0] == ' ')
        s.erase(s.begin());
    while (s.back() == ' ')
        s.erase(s.end() - 1);
}

vector<string> tokenize(string &s)
{
    vector<string> res;
    string word = "";
    int l = s.length();
    for (int j = 0; j <= l; ++j)
    {
        if (j == l or s[j] == ' ')
        {
            if (word != "")
                res.push_back(word);
            word = "";
        }
        else
            word = word + s[j];
    }
    return res;
}