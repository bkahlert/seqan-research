#include <iostream>
#include <seqan/index.h>

using namespace seqan;

template <typename TIndexIt>
void DFS(TIndexIt it){
	std::cout << representative(it) << std::endl;
	if (goDown(it))
		do
		DFS(it);
	while (goRight(it));

}
int main()
{
	typedef Index<CharString> TIndex;
	TIndex index("tobeornottobe");
	Iterator< TIndex, TopDown<> >::Type it(index);
	DFS(it);
	return 0;
}