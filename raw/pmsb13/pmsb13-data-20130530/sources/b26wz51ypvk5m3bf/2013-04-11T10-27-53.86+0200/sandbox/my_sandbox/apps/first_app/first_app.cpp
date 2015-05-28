#include <seqan/sequence.h>
#include <seqan/index.h>
#include <seqan/basic.h>
int main ()
{
    typedef Index<DnaString, IndexQGram< UngappedShape<3> > > TIndex;
    TIndex index("CATGATTACATA");
    hash(indexShape(index), "CAT");
    for (unsigned i = 0; i < seqan::length(getOccurrences(index, indexShape(index))); ++i)
        std::cout << seqan::getOccurrences(index, indexShape(index))[i] << std::endl;

    return 0;
}
