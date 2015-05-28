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

	Finder<CharString> finder(haystack);
    Pattern<StringSet<CharString>, SetHorspool> pattern(needles);

    while (find(finder, pattern))
        std::cout << '[' << beginPosition(finder) << ',' << endPosition(finder) << ")\t" << infix(finder) << std::endl;
    return 0;
}
