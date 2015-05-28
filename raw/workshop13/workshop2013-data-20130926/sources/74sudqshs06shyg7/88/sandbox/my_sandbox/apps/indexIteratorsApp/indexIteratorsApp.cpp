#include <iostream>
#include <seqan/index.h>

using namespace seqan;
int main ()
{
    typedef Index<CharString> TIndex;  
    TIndex index("mariaM");
    Iterator< TIndex, TopDown<Preorder> >::Type it(index);
    
    std::cout << representative(it) << std::endl;
    while(!atEnd(it))
        std::cout << representative(it) << std::endl;
        ++it;        
    } 
    return 0;
}