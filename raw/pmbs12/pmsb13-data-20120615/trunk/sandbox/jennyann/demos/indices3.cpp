#include <iostream>
#include <seqan/index.h>

using namespace seqan;

int main ()
{
    String<char> myString = "tobeornottobe";
    typedef Index< String<char> > TMyIndex;
    TMyIndex myIndex(myString);

	Iterator< TMyIndex, TopDown< ParentLinks<Preorder> > >::Type myIterator(myIndex);
    
	do {
		std::cout << representative(myIterator) << std::endl;
		if (!goDown(myIterator) && !goRight(myIterator))
			 while (goUp(myIterator) && !goRight(myIterator)) ;
	} while (!isRoot(myIterator));
    return 0;
}