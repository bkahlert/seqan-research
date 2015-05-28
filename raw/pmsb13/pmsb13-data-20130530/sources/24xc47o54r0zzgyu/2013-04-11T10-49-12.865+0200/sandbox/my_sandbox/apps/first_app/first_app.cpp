#include <seqan/sequence.h>
#include <seqan/index.h>

using namespace seqan;

int main()
{
	typedef Index<CharString> TIndex;
	TIndex index("tobeornottobe");
	Iterator< TIndex, TopDown<ParentLinks<> > >::Type it(index);

	do
	{
		std::cout<<representative(it)<<std::endl;
		while(goDown(it));
		while(!goRight(it))
			goUp(it);

	}
	while(!isRoot(it));

    return 0;
}