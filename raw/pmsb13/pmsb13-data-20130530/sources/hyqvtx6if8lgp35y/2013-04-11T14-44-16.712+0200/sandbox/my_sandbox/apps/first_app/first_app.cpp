#include <iostream>
#include <seqan/index.h>
#include <seqan/find.h>

using namespace seqan;


int main()
{
	CharString haystack = "tobeornottobe";
    CharString pattern = "tobe";
	Index<CharString, IndexQGram<UngappedShape<4>, OpenAddressing> > index(haystack);
    Finder<CharString> finder(index);
    while (find(finder, pattern))
	{	            
		std::cout << '[' << beginPosition(finder) << ',';
 		std::cout << endPosition(finder) << ")\t" << infix(finder) << std::endl;
	}
	return 0;
}
