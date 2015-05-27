#include <iostream>
#include <seqan/index.h>

using namespace seqan;


template < typename TIndexSpec >
void indexfinder()
{
    String<char> myString = "tobeornottobe";
    typedef Index< String<char>, TIndexSpec > TMyIndex;
    TMyIndex myIndex(myString);

	Iterator< TMyIndex, TopDown< ParentLinks<Preorder> > >::Type myIterator(myIndex);
    
	do {
		std::cout << representative(myIterator) << std::endl;
		do {
		if (!goDown(myIterator) && !goRight(myIterator))
        while (goUp(myIterator) && !goRight(myIterator)) ;
		} while (repLength(myIterator) > 3);
	} while (!isRoot(myIterator));
	std::cout << std::endl;
}

int main ()
{
	indexfinder< IndexEsa<> > ();
	indexfinder< IndexWotd<> > ();
    return 0;
}