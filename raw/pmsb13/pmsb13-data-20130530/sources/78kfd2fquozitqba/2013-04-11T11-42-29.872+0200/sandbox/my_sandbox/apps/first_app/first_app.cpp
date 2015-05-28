#include <iostream>
#include <seqan/index.h>

using namespace seqan;
int main ()
{
     String<Dna> text = "CATGATTACATA";
     Index<DnaString, IndexQGram< OneGappedShape > > index(text);
     stringToShape(indexShape(index),"1101");
     hash(indexShape(index), "AT-A");
     for (unsigned i = 0; i < length(getOccurrences(index, indexShape(index))); ++i)
         std::cout << getOccurrences(index, indexShape(index))[i] << std::endl;

     return 0;
 }

