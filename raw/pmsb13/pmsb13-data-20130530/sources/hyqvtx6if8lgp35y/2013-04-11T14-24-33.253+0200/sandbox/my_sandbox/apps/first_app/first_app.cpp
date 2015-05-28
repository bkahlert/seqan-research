#include <iostream>
#include <seqan/index.h>
#include <seqan/find.h>

using namespace seqan;


int main()
{
	CharString haystack = "Simon, send more money!";
    CharString needle = "more";

    Finder<CharString> finder(haystack);
    Pattern<CharString, Myers< > > pattern(needle);
    while (find(finder, pattern, -2))
	{
        while (findBegin(finder, pattern, getScore(pattern)))
		{		            
			std::cout << '[' << beginPosition(finder) << ',';
 			std::cout << endPosition(finder) << ")\t" << infix(finder) << std::endl;
		}
	}
	return 0;
}
