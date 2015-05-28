#include <seqan/sequence.h>
#include <seqan/index.h>

using namespace seqan;

int main()
{
	typedef Index<CharString> TIndex;
	TIndex index("tobeornottobe");
	Iterator<TIndex,TopDown<>>::Type it(index);

	std::cout << representative(it) << std::endl;
	goDown(it);
	std::cout << representative(it) << std::endl;

	return 0;

}