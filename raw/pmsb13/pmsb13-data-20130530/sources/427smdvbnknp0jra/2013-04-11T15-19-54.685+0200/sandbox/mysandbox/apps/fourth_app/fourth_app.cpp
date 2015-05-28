#include <seqan/sequence.h>
#include <seqan/index.h>

using namespace seqan;

template <typename Indexinput>
void showTree(String<char> & input)
{	
	std::cout << std::endl;

	typedef Index<String<char>,Indexinput> TIndex;
	TIndex index(input);	
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
	Iterator<StringSet<String<char>>>::Type input;
	appendValue(input,"tobeornottobe");
	appendValue(input,"thebeeonthecomb");
	appendValue(input,"beingjohnmalkovich");

	showTree<IndexEsa<>>(*input);
	
	showTree<IndexWotd<>>("tobeoornottobe");
	
	
	return 0;

}