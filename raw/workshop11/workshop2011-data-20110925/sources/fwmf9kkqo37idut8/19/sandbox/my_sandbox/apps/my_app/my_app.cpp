#include <iostream>
#include <seqan/sequence.h>  // CharString, ...
#include <seqan/file.h>      // to stream a CharString into cout
#include <seqan/basic.h>
using namespace seqan;

void tut01();

int main(int, char **) {
	tut01();

	return 0;
}

template <typename TAlphabet>
void showAllLetterOfMyAlphabet(TAlphabet const &)
{
	typedef typename Size<TAlphabet>::Type TSize;
	TSize alphSize = ValueSize<TAlphabet>::VALUE;
	for (TSize i = 0; i < alphSize; ++i)
		std::cout << i << ',' << TAlphabet(i) << "  ";
	std::cout << std::endl;
}

void tut01(){
	showAllLetterOfMyAlphabet(AminoAcid());
	showAllLetterOfMyAlphabet(Dna());
	showAllLetterOfMyAlphabet(Dna5());
}