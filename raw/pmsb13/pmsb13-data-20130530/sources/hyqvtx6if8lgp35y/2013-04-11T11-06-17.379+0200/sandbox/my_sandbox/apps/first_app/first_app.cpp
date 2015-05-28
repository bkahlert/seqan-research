#include <iostream>
#include <seqan/index.h>

using namespace seqan;

int main ()
{
    typedef Index<CharString> TIndex;
    TIndex index("tobeornottobe");
	int count = 0;
    Iterator< TIndex, TopDown<ParentLinks<> > >::Type it(index);    
	do 
	{
		std::cout << count << '\t' << representative(it) << std::endl;	
  		++count;
		if (!goDown(it) && !goRight(it))
			while (goUp(it) && !goRight);
    } while (!isRoot(it));
    
	return 0;
}
