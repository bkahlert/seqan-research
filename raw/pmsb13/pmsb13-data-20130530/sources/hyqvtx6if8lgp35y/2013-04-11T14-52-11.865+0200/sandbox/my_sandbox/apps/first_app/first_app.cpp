#include <iostream>
#include <seqan/index.h>
#include <seqan/find.h>

using namespace seqan;

typedef Index<CharString, IndexQGram<UngappedShape<4>, OpenAddressing> > TIndex;

int main()
{
	CharString haystack = "tobeornottobe";
    CharString pattern = "tobe";
	TIndex index(haystack);
    Finder<TIndex> finder(index);
	
    while (find(finder, pattern))
	{	            
		std::cout << '[' << beginPosition(finder) << ',';
 		std::cout << endPosition(finder) << ")\t" << infix(finder) << std::endl;
		
	}
	clear(finder);
	return 0;
}
