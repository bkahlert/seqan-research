// FRAGMENT(includes)
#include <iostream>
#include <seqan/index.h>

using namespace seqan;
using namespace std;
int main ()
{
    Index<DnaString, IndexQGram<OneGappedShape> > index("CATGATTACATA");
    stringToShape(indexShape(index), "1101");
    hash(indexShape(index), "AT-A");
    for (unsigned i = 0; i < length(getOccurrences(index, indexShape(index))); ++i)
        cout << getOccurrences(index, indexShape(index))[i] << endl;

    return 0;
}
