#include <seqan/sequence.h>
#include <seqan/index.h>

using namespace seqan;

char[] result;

template <typename Indexinput,typename T>
void showTree(T input)
{	
	clear(result);
	
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
		appendValue(result,";");
	}
}

void lookOf()
{	
	StringSet<char> t;
	String<Repeat<unsigned, unsigned> > repeats;

	findRepeats(repeats, result, 3);
	//std::cout << repeats;
	/*
	for(unsigned i=0; i<length(result);i++)
	{
		std::cout << result[i] <<std::endl;	
	}*/

}

int main()
{
	StringSet<String<char>> input;
	appendValue(input,"tobeornottobe");
	appendValue(input,"thebeeonthecomb");
	appendValue(input,"beingjohnmalkovich");

	showTree<IndexEsa<>>(input);
	lookOf();
	//showTree<IndexWotd<>>("tobeoornottobe");
	
	
	return 0;

}