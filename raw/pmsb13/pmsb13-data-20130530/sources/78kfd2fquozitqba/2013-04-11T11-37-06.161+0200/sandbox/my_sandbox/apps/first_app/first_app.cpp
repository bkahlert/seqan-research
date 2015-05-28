#include <iostream>
#include <seqan/index.h>

using namespace seqan;
int main ()
{
     Index<DnaString, IndexQGram< OneGappedShape > > index;
     index("CATGATTACATA");
     hash(indexShape(index), "AT-A");
    return 0;
}

