#include <iostream>
#include <seqan/sequence.h>  // CharString, ...
#include <seqan/file.h>      // to stream a CharString into cout
#include <seqan/basic.h>
using namespace seqan;

int main(int, char **) {
	std::cout << "Hello World!" << std::endl;
	seqan::CharString mySeqanString = "Hello SeqAn!";
	std::cout << mySeqanString << std::endl;
	return 1;
}

void tut01(){
	template <typename TAlphabet>
	void showAllLetterOfMyAlphabet(TAlphabet const &)
	{
		typedef typename Size<TAlphabet>::Type TSize;
		TSize alphSize = ValueSize<TAlphabet>::VALUE;
		for (TSize i = 0; i < alphSize; ++i)
			std::cout << i << ',' << TAlphabet(i) << "  ";
		std::cout << std::endl;
	}

	showAllLetterOfMyAlphabet(AminoAcid());
	showAllLetterOfMyAlphabet(Dna());
	showAllLetterOfMyAlphabet(Dna5());
	return 0;

}