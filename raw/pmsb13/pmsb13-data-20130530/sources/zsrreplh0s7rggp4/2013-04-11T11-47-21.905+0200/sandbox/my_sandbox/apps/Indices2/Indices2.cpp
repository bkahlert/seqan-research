#include <iostream>
#include <seqan/index.h>

using namespace seqan;
int main ()
{
    Index<DnaString, IndexQGram<OneGappedShape> > index("CATGATTACATA");
//  Index<DnaString, IndexQGram<GenericShape> > index("CATGATTACATA");
    stringToShape(indexShape(index), "1101");
//  Index<DnaString, IndexQGram<GappedShape<HardwiredShape<1,2> > > > index("CATGATTACATA");
     hash(indexShape(index), "ATCA");
    for (unsigned i = 0; i < length(getOccurrences(index, indexShape(index))); ++i)
        std::cout << getOccurrences(index, indexShape(index))[i] << std::endl;

    return 0;
}
