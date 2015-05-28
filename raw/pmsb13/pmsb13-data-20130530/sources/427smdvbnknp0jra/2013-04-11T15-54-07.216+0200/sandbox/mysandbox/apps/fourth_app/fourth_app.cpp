#include <seqan/sequence.h>
#include <seqan/index.h>

using namespace seqan;

StringSet<String<char>> result;

template <typename Indexinput,typename T>
void showTree(T input)
{	
	StringSet<String<char>> result;
	
	Iterator<T>::Type it1 = begin(input);
	
	for(;!atEnd(it1);it1++)
	{
		std::cout << std::endl;
	
		typedef Index<String<char>,Indexinput> TIndex;
		TIndex index(*it1);	
		Iterator<TIndex,TopDown<ParentLinks<>>>::Type it2(index);
		
		appendValue(result,representative(it2));
		goDown(it2);
		
		while(!isRoot(it2))
		{	
	
			appendValue(result,representative(it2));
				
			if (!goDown(it2) && !goRight(it2))
			{
				while (goUp(it2) && !goRight(it2)){}	
			}
		}
	}
}

void lookOf()
{


}

int main()
{
	StringSet<String<char>> input;
	appendValue(input,"tobeornottobe");
	appendValue(input,"thebeeonthecomb");
	appendValue(input,"beingjohnmalkovich");

	showTree<IndexEsa<>>(input);
	
	//showTree<IndexWotd<>>("tobeoornottobe");
	
	
	return 0;

}