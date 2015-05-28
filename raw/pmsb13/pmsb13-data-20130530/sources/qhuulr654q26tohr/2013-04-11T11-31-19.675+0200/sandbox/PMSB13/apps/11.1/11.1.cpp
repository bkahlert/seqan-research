#include <seqan/sequence.h>
#include <seqan/index.h>

using namespace seqan;
using namespace std;

int main()
{
    typedef String<Dna> TAlphabet;
    TAlphabet haystack = "CATGATTACATA";
    TAlphabet needle = "AT-A";
    
    typedef Index<TAlphabet, IndexQGram< OneGappedShape<> > > TIndex;
    TIndex index(haystack);
    hash(indexShape(index)("1101"), "CAT");
    
    return 0;
}