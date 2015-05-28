#include <iostream>
#include <seqan/find.h>
#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/modifier.h>

using namespace seqan;

int main()
{
    String<char> haystack = "Simon, send more money!";
    String<String<char>> needle;
	appendValue(needle,"more");
	appendValue(needle,"send");
	appendValue(needle,"mo");
    Finder<String<char>> finder(haystack);
    Pattern<String<String<char>>, WuManber> pattern(needle);
    while (find(finder, pattern))
		std::cout << '[' << beginPosition(finder) << ',' << endPosition(finder) << ")\t" << infix(finder) << std::endl;
	return 0;
}