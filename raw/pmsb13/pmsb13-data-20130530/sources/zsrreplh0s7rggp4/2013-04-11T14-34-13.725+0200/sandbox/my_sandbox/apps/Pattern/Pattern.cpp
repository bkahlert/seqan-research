#include <iostream>
#include <seqan/find.h>

using namespace seqan;

int main()
{
    CharString haystack = "Simon, send more money!";
    CharString needle="more";
    
    Finder<CharString> finder(haystack);
    Pattern<CharString, Myers> pattern(needle,SimpleScore(0, -2, -1));
    while (find(finder, pattern, -2)){
        std::cout << '[' << beginPosition(finder) << ',' << endPosition(finder) << ")\t" << infix(finder) << "\t";
        std::cout << position(pattern)  << '\t' <<  infix(finder) << std::endl;
    }
    return 0;
}