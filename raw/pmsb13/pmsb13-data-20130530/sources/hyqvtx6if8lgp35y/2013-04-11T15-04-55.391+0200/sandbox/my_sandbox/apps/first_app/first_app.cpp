#include <iostream>
#include <seqan/index.h>
#include <seqan/find.h>

using namespace seqan;

//typedef Index<CharString, IndexQGram<UngappedShape<4>, OpenAddressing> > TIndex;
typedef Index<CharString > TIndex;
int main()
{
	CharString haystack = "tobeornottobe";
    CharString pattern = "tobe";
	TIndex index(haystack);
    Finder<TIndex> finder(index);
	//clear(finder)

	int i = 0;
    while (find(finder, pattern))
	{	            
	
		std::cout << '[' << beginPosition(finder) << ',';
 		std::cout << endPosition(finder) << ")\t" << infix(finder) << std::endl;
		std::cout << i << std::endl;
		++i;
	}
	
	return 0;
}
