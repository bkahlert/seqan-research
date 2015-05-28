#include <iostream>
#include <seqan/index.h>

using namespace seqan;
int main ()
{
     Index<DnaString, IndexQGram< OneGappedShape<3> > > TIndex = index("CATGATTACATA");
    return 0;
}

