#include <iostream>
#include <seqan/find.h>
#include <seqan/index.h>

using namespace seqan;

int main()
{
    Index<CharString> index ("tobeornottobe");
    CharString needle = "be";
    Finder<Index<CharString> > finder (index);

	Pattern<CharString> pattern (needle);
    while (find(finder, pattern))
        std::cout << '[' << beginPosition(finder) << ',' << endPosition(finder) << ")\t" << infix(finder) << std::endl;

	clear(finder);
    while (find(finder, "be"))
        std::cout << '[' << beginPosition(finder) << ',' << endPosition(finder) << ")\t" << infix(finder) << std::endl;

    return 0;
}