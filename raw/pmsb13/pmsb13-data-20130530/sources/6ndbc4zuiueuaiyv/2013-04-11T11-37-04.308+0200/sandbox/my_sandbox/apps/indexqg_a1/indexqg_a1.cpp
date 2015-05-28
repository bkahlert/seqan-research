#include <iostream>
#include <seqan/index.h>

using namespace seqan;
using namespace std;

int main ()
{
    typedef Index<DnaString, IndexQGram< OneGappedShape<3> > > TIndex;
    TIndex index("CATGATTACATA");
    
    hash(indexShape(index), "ATA");
    for (unsigned i = 0; i < length(getOccurrences(index, indexShape(index))); ++i)
        std::cout << getOccurrences(index, indexShape(index))[i] << std::endl;

    return 0;
}

