#include <iostream>
#include <seqan/find.h>

using namespace seqan;
using namespace std;

int main()
{
    CharString haystack = "Simon, send more money!";
    StringSet< String<char> > needles;
    appendValue(needles, "mo");
    appendValue(needles, "send");
    appendValue(needles, "more");
    
    Finder<CharString> finder(haystack);

    cout << "Pattern: " << needleSet[i] << endl;
    Pattern<CharString, Horspool> pattern(needleSet[i]);
    while (find(finder, pattern))
        std::cout << '[' << beginPosition(finder) << ',' << endPosition(finder) << ")\t" << infix(finder) << std::endl;
    clear(finder);

    return 0;
}
