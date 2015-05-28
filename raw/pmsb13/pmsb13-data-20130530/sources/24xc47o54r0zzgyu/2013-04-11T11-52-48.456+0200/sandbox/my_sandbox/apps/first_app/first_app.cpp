#include <seqan/sequence.h>
#include <seqan/index.h>

using namespace seqan;

int main()
{
	Shape<CharString> shape;
	stringToShape(shape, "1101");
	typedef Index<CharString, IndexQGram< GappedShape<HardwiredShape<0,1> > > > TIndex;
    TIndex index("CATGATTACATA");
	hash(shape, "AT-A");

    for (unsigned i = 0; i < length(getOccurrences(index, indexShape(index))); ++i)
        std::cout << getOccurrences(index, indexShape(index))[i] << std::endl;


    return 0;
}