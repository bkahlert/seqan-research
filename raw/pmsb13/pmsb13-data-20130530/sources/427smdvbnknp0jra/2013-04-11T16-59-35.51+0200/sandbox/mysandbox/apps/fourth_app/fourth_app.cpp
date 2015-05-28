#include <seqan/sequence.h>
#include <seqan/index.h>

using namespace seqan;

template <typename Indexinput,typename T>
void showTree(T input)
{	
	typedef Index<String<char>,Indexinput> TIndex;
	TIndex index(input);	
	Iterator<TIndex,TopDown<ParentLinks<Preorder>>>::Type it2(index);
	
	
	while(!atEnd(it2))
	{	
		std::cout << representative(it2) << std::endl;
		it2++;
	}
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