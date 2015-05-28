#include <iostream>
#include <seqan/find.h>

using namespace seqan;

int main()
{
    CharString haystack = "Simon, send more money!";
    StringSet< String<char> > needleSet;
    resize(needleSet,3);
    needleSet[0] = "mo";
    needleSet[1] = "send";
    needleSet[2] = "more";
    
    Finder<CharString> finder(haystack);
    for (unsigned i = 0; i < length(needleSet); ++i)
    {
        Pattern<CharString, Horspool> pattern(needleSet[i]);
        while (find(finder, pattern))
            std::cout << '[' << beginPosition(finder) << ',' << endPosition(finder) << ")\t" << infix(finder) << std::endl;
    }
    return 0;
}
