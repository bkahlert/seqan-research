#include <seqan/sequence.h>
#include <seqan/index.h>

using namespace seqan;

int main()
{
	typedef Index<CharString> TIndex;
	TIndex index("tobeornottobe");
	Iterator<TIndex,TopDown<ParentLinks<>>>::Type it(index);
	
	std::cout << representative(it) << std::endl;
	goDown(it);
	
	while(!isRoot(it))
	{
	    std::cout << representative(it) << std::endl;
		if (!goDown(it) && !goRight(it))
		{
            while (goUp(it) && !goRight(it)) ;
		}
	}
	
	return 0;

}