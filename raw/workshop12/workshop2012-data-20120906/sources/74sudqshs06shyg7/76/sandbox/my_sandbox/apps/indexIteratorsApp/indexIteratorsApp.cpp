#include <iostream>
#include <seqan/index.h>

using namespace seqan;
int main ()
{
    typedef Index<CharString> TIndex;  
    TIndex index("mariaM");
    Iterator< TIndex, TopDown<PreOrder > >::Type it(index);
    
    
    do {
        std::cout << representative(it) << std::endl;
        if (!goDown(it) && !goRight(it))
            while (goUp(it) && !goRight(it)) ;
            
    } while (!isRoot(it));
    return 0;
}