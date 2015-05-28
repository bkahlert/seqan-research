#include <iostream>
#include <seqan/find.h>

using namespace seqan;

int main()
{
    CharString haystack = "Simon, send more money!";
	StringSet<CharString> needles;
	appendValue(needles,  "mo");
	appendValue(needles,  "send");
	appendValue(needles,  "more");

    //CharString needle = "mo";
	Finder<CharString> finder(haystack);
    Pattern<String<CharString>, SetHorspool> pattern(needles);
	//setNeedle(pattern,needles)
    while (find(finder, pattern))
        std::cout << '[' << beginPosition(finder) << ',' << endPosition(finder) << ")\t" << infix(finder) << std::endl;
    return 0;
}
