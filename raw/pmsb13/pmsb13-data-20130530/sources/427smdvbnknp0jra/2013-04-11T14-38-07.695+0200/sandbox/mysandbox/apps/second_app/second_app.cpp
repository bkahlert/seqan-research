#include <seqan/sequence.h>
#include <seqan/index.h>

using namespace seqan;

void showTree(int input)
{	
	if(input==1)
	{
		typedef Index<String<char>,IndexEsa<>> TIndex;
	}
	else
	{
		typedef Index<String<char>,IndexWotd<>> TIndex;
	}
	TIndex index("tobeornottobe");
	Iterator<TIndex,TopDown<ParentLinks<>>>::Type it(index);
	
	std::cout << representative(it) << std::endl;
	goDown(it);
	
	while(!isRoot(it))
	{	

		if(repLength(it)<=3)
		{
			std::cout << representative(it) << std::endl;
		}

		if (!goDown(it) && !goRight(it))
		{
			while (goUp(it) && !goRight(it)){}
				
		}
	}
}

int main()
{

	showTree(IndexEsa<>);

}