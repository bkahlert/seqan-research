#include <iostream>
#include <seqan/index.h>

using namespace seqan;
int main ()
{
    typedef Index<CharString> TIndex;
    TIndex index("tobeornottobe");
    Iterator< TIndex, TopDown<ParentLinks<> > >::Type it(index);
	unsigned i = 0;
    do {
        std::cout << representative(it) << std::endl;
		++i;
        if (!goDown(it) && !goRight(it))
		{
			if(!goRight(it))
				while (goUp(it) && !goRight(it)) 
					--i;
		}
		else ++i;
           
    } while (!isRoot(it));
    return 0;
}