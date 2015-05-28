#include <seqan/sequence.h>
#include <seqan/index.h>

using namespace seqan;

int main()
{
	typedef Index<DnaString, IndexQGram< GappedShape<HardwiredShape<> > > > TIndex;
    TIndex index("CATGATTACATA");

	hash(indexShape(index), "AT-A");
    for (unsigned i = 0; i < length(getOccurrences(index, indexShape(index))); ++i)
        std::cout << getOccurrences(index, indexShape(index))[i] << std::endl;


    return 0;
}