#include <iostream>
#include <seqan/find.h>

using namespace seqan;

int main()
{
    CharString haystack = "Simon, send more money!";
    StringSet<CharString> needle;
	appendValue(needle, "mo");
	appendValue(needle, "send");
	appendValue(needle, "more");

	Finder<CharString> finder(haystack);
    Pattern<CharString, WuManber> pattern(needle);
    while (find(finder, pattern))
        std::cout << '[' << beginPosition(finder) << ',' << endPosition(finder) << ")\t" << infix(finder) << std::endl;

	return 0;

}