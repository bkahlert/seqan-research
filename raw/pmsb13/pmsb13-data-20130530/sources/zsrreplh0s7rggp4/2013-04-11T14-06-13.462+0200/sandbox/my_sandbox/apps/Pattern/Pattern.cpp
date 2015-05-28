#include <iostream>
#include <seqan/find.h>

using namespace seqan;

int main()
{
    CharString haystack = "Simon, send more money!";
    String<CharString> needle;
    
    appendValue(needle,"mo");
    appendValue(needle,"money");
    appendValue(needle,"n");
    
    Finder<CharString> finder(haystack);
    Pattern<String<CharString>, WuManber> pattern(needle);
    while (find(finder, pattern))
        std::cout << '[' << beginPosition(finder) << ',' << endPosition(finder) << ")\t" << infix(finder) << std::endl;
        std::cout << position(pattern)  << '\t' <<  infix(finder) << std::endl;
    return 0;
}