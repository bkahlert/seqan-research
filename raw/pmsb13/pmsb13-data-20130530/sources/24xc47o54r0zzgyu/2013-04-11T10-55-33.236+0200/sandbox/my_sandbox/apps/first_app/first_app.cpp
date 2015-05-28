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
		while(goDown(it)) std::cout<<"down"<<std::endl;	
		while(!goRight(it))
		{
			goUp(it);
			std::cout<<"right"<<std::endl;
		}

	} while(!isRoot(it));

    return 0;
}