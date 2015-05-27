#include <iostream>
#include <seqan/index.h>

using namespace seqan;

int main ()
{
	typedef String<char> TSequence;
	TSequence myString1 = "CDFGHC";
	TSequence myString2 = "CDEFGAHC";

	StringSet<TSequence> sequences;

	appendValue(sequences,myString1);
    appendValue(sequences,myString2);

	typedef Index< StringSet<TSequence> > TMyIndex;
    TMyIndex myIndex(sequences);

	 Iterator< TMyIndex, Mums >::Type myIterator(myIndex);
    //goBegin(myIterator);
    while (!atEnd(myIterator))
    {
        std::cout << representative(myIterator) << std::endl;
        ++myIterator;
    }
    return 0;
}