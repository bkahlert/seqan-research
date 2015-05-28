#include <iostream>
#include <seqan/index.h>

using namespace seqan;
using namespace std;

void nodeTraverser(indexType)
{
    typedef Index<CharString> TIndex;
    TIndex index("tobeornottobe");
    typename Iterator< TIndex, TopDown<ParentLinks<> > >::Type it(index);


    do {
        std::cout << representative(it) << std::endl;
        if (!goDown(it) || repLength(it) > 3)
            do {
                if (length(representative(it)) <= 3)
                    cout << representative(it) << endl;
                if (!goDown(it) && !goRight(it))
                } while (goUp(it) && !goRight(it));
    } while (!isRoot(it));
    std::cout << std::endl;
}


int main ()
{
    nodeTraverser(
    return 0;
}

