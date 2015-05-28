#include <iostream>
#include <seqan/index.h>

using namespace seqan;
int main ()
{
     String<Dna> text = "CATGATTACATA";
     Index<DnaString, IndexQGram< OneGappedShape > > index(text);
     hash(indexShape(index), "AT-A");
    return 0;
}

