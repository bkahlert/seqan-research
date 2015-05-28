#include <seqan/sequence.h>
#include <seqan/index.h>
#include <iostream>

using namespace std;
using namespace seqan;

int main() {
	String < Dna > text = "TTATTAAGCGTATAGCCCTATAAATATAA";
	String < Dna > pattern = "TATAA";
	Index < String<Dna> , IndexEsa<> > esaIndex(text);
	Finder < Index<String<Dna> , IndexEsa<> > > esaFinder(esaIndex);

	while (find(esaFinder, pattern)) {
		cout << "Found " << pattern << " at position " << position(esaFinder)
				<< endl;
	}

	CharString text2 = "tobeornottobe";
	Index < CharString, IndexEsa<> > esaindex2(text2);
	typedef Iterator<Index<CharString, IndexEsa<> > , TopDown<ParentLinks<
			Preorder> > >::Type TIterator;

	TIterator it(esaindex2);
	do {
		cout << representative(it) << endl;
		if (!goDown(it) && !goRight(it))
			while (goUp(it) && !goRight(it))
				;
	} while (!root(it));

	return 0;
}
