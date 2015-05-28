#include <seqan/sequence.h>
#include <seqan/index.h>

using namespace seqan;

template <typename Indexinput,typename T>
void showTree(T input)
{	
	Iterator<T>::Type it1 = begin(input);
	std::cout << it1 <<std::endl;
	for(;!atEnd(it1);i++)
	{
		std::cout << it1 <<std::endl;
		std::cout << std::endl;
	
		typedef Index<String<char>,Indexinput> TIndex;
		TIndex index(*it1);	
		Iterator<TIndex,TopDown<ParentLinks<>>>::Type it2(index);
		
		std::cout << representative(it2) << std::endl;
		goDown(it2);
		
		while(!isRoot(it2))
		{	
	
			if(repLength(it2)<=3)
			{
				std::cout << representative(it2) << std::endl;
			}
	
			if (!goDown(it2) && !goRight(it2))
			{
				while (goUp(it2) && !goRight(it2)){}
					
			}
		}

	}
	
}

int main()
{
	StringSet<String<char>> input;
	appendValue(input,"tobeornottobe");
	appendValue(input,"thebeeonthecomb");
	appendValue(input,"beingjohnmalkovich");

	showTree<IndexEsa<>>(input);
	
	showTree<IndexWotd<>>("tobeoornottobe");
	
	
	return 0;

}