#include <iostream>
#include <seqan/index.h>

using namespace seqan;

int main ()
{
	typedef String<char> TSequence;
	TSequence myString1 = "tobeornottobe";
	TSequence myString2 = "thebeeonthecomb";
	TSequence myString3 = "beingjohnmalkovich";

	StringSet<TSequence> sequences;

	appendValue(sequences,myString1);
    appendValue(sequences,myString2);
	appendValue(sequences,myString3);

	typedef Index< StringSet<TSequence> > TMyIndex;
    TMyIndex myIndex(sequences);

	 Iterator< TMyIndex, TopDown< ParentLinks<Postorder> > >::Type myIterator(myIndex);
    goBegin(myIterator);
    while (!atEnd(myIterator))
    {
        std::cout << representative(myIterator) << std::endl;
        ++myIterator;
    }
    return 0;
}