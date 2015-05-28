#include <iostream>
#include <seqan/index.h>

using namespace seqan;
int main ()
{
    typedef Index<CharString> TIndex;  
    TIndex index("mariaM");
    Iterator< TIndex, TopDown<Preorder> >::Type it(index);
    
    std::cout << representative(it) << std::endl;
 /*  do {
        std::cout << representative(it) << std::endl;
        ++it;
      /*  if (!goDown(it) && !goRight(it))
            while (goUp(it) && !goRight(it)) ;
            
    } while (!isRoot(it));*/
    return 0;
}