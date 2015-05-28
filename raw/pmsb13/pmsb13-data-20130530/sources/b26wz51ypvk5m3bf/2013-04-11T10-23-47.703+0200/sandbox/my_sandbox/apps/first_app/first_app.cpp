#include <seqan/sequence.h>
#include <seqan/index.h>
int main ()
{
    typedef Index<DnaString, IndexQGram< UngappedShape<3> > > TIndex;
    TIndex index("CATGATTACATA");
    hash(indexShape(index), "CAT");
    for (unsigned i = 0; i < length(getOccurrences(index, indexShape(index))); ++i)
        std::cout << seqan::getOccurrences(index, indexShape(index))[i] << std::endl;

    return 0;
}
