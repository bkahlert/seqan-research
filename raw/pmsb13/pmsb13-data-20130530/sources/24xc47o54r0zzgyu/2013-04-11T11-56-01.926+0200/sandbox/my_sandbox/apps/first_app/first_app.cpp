#include <seqan/sequence.h>
#include <seqan/index.h>

using namespace seqan;

int main()
{
	Shape<DnaString> shape;
	typedef Index<DnaString, IndexQGram< OneGappedShape > > TIndex;
    TIndex index("CATGATTACATA");
	stringToShape(shape, "1101");
	hash(shape, "ATCA");

    for (unsigned i = 0; i < length(getOccurrences(index, shape)); ++i)
        std::cout << getOccurrences(index, shape)[i] << std::endl;


    return 0;
}