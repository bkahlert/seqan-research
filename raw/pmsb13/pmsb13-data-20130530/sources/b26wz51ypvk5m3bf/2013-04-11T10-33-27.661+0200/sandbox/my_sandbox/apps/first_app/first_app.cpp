#include <seqan/sequence.h>
#include <seqan/index.h>
#include <seqan/basic.h>
#include <iostream>
#include <seqan/file.h>
#include <seqan/modifier.h>
int main ()
{
    typedef Index<DnaString, IndexQGram< UngappedShape<3> > > TIndex;
    TIndex index("CATGATTACATA");
    seqan::hash(indexShape(index), "CAT");
    for (unsigned i = 0; i < length(seqan::getOccurrences(index, indexShape(index))); ++i)
        std::cout << getOccurrences(index, indexShape(index))[i] << std::endl;

    return 0;
}
