#include <iostream>
#include <seqan/find.h>
#include <seqan/index.h>

using namespace seqan;

int main()
{
	typedef StringSet<CharString> THaystacks; 
	THaystacks haystacks;
	appendValue(haystacks, "tobeornottobe");
	appendValue(haystacks, "thebeeonthecomb");
	appendValue(haystacks, "beingjohnmalkovich");

    Index<THaystacks> index (haystacks);
    CharString needle = "be";
    Finder<Index<THaystacks> > finder (index);

	Pattern<CharString> pattern (needle);
    while (find(finder, pattern))
        std::cout << '[' << beginPosition(finder) << ',' << endPosition(finder) << ")\t" << infix(finder) << std::endl;

	clear(finder);
    while (find(finder, "be"))
        std::cout << '[' << beginPosition(finder) << ',' << endPosition(finder) << ")\t" << infix(finder) << std::endl;

    return 0;
}