#include <seqan/sequence.h>
#include <seqan/index.h>

using namespace seqan;
using namespace std;

int main()
{
    typedef String<Dna> TAlphabet;
    TAlphabet haystack = "CATGATTACATA";
    TAlphabet needle = "AT-A";
    
    typedef Index<TAlphabet, IndexQGram< OneGappedShape> > TIndex;
    TIndex index(haystack);
    stringToShape(indexShape(index),"1101");
    Iterator<TAlphabet, Rooted>::Type it = begin(needle);
    hash(indexShape(index), it);
    
    for (unsigned i = 0; i < length(getOccurrences(index, indexShape(index))); ++i)
	std::cout << getOccurrences(index, indexShape(index))[i] << std::endl;
    
    return 0;
}