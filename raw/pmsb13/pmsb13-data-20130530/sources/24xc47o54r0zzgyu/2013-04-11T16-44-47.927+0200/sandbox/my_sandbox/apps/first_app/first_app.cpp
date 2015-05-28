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
		if (i<=3)
		{
			if (!goDown(it))
			{
				if(!goRight(it))
				{
					while (!goRight(it))
					{
						goUp(it);
						--i;
					}
				}
			}
			else ++i;
		}
		else 
		{
			while (!goRight(it))
					{
						goUp(it);
						--i;
					}
		}
    } while (!isRoot(it));
    return 0;
}