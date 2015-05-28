#include <iostream>
#include <seqan/find.h>

using namespace seqan;

int main()
{
    CharString haystack = "Simon, send more money!";
    String<CharString> needle;
	appendValue(needle,"mo");
	appendValue(needle,"more");
	appendValue(needle,"send");
    Finder<CharString> finder(haystack);
    Pattern<CharString, WoManber> pattern(needle);
    while (find(finder, pattern))
        std::cout << '[' << beginPosition(finder) << ',' << endPosition(finder) << ")\t" << infix(finder) << std::endl;

    return 0;
}

