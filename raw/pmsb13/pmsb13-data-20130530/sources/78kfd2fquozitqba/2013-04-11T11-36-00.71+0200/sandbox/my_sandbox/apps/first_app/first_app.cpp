#include <iostream>
#include <seqan/index.h>

using namespace seqan;
int main ()
{
     Index<DnaString, IndexQGram< OneGappedShape > > index;
    return 0;
}

