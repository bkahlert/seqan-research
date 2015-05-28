#include <iostream>
#include <seqan/find.h>

using namespace seqan;

int main()
{
	// Assignment 1
    CharString haystack = "Simon, send more money!";
	StringSet<CharString> needles;
	appendValue(needles,  "mo");
	appendValue(needles,  "send");
	appendValue(needles,  "more");

	Finder<CharString> finder(haystack);
    Pattern<StringSet<CharString>, MultipleShiftAnd> pattern(needles);

    while (find(finder, pattern))
        std::cout << '[' << beginPosition(finder) << ',' << endPosition(finder) << ")\t" << position(pattern) << std::endl;
	
	std::cout << "Assignment 2" << std::endl; 

	// Assignment 2
	CharString needle = "more";
	Pattern<CharString, Myers<CharString, EditDistance> > pattern2(needle);
    while (find(finder, pattern2/*, getScore(pattern2)*/))
       // while (findBegin(finder, pattern2, getScore(pattern2)))
            std::cout << '[' << beginPosition(finder) << ',' /*<< endPosition(finder)*/ << ")\t" << infix(finder) << std::endl;



    return 0;
}
