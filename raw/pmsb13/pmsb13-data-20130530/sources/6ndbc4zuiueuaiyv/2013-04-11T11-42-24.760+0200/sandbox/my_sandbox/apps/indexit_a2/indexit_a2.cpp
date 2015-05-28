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
        if (!goDown(it) && !goRight(it))
            while (goUp(it) && !goRight(it));
    } while (!isRoot(it));
    return 0;
}

