#include <iostream>
#include <seqan/index.h>

using namespace seqan;
using namespace std;

int main ()
{
    typedef Index<CharString> TIndex;
    TIndex index("tobeornottobe");
    Iterator< TIndex, TopDown<ParentLinks<> > >::Type it(index);
    do {
        cout << representative(it) << endl;
    } while (!isRoot(it));
    return 0;
}

