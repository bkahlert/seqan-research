#include <seqan/sequence.h>
#include <seqan/index.h>
#include <seqan/basic.h>
#include <iostream>
#include <seqan/file.h>
#include <seqan/modifier.h>
int main ()
{
    typedef seqan::Index<seqan::DnaString, seqan::IndexQGram< seqan::GappedShape<4> > > TIndex;
    TIndex index("CATGATTACATA");
    seqan::hash(indexShape(index), "AT-A");
    for (unsigned i = 0; i < seqan::length(seqan::getOccurrences(index, indexShape(index))); ++i)
        std::cout << seqan::getOccurrences(index, indexShape(index))[i] << std::endl;

    return 0;
}
