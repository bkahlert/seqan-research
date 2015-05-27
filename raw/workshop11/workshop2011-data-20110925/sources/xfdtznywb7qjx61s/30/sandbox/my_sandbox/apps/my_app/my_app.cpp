#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <iostream>
using namespace seqan;

template <typename TAlphabet>
void showAllLetterOfMyAlphabet(TAlphabet const &)
{
    typedef typename Size<TAlphabet>::Type TSize;
    TSize alphSize = ValueSize<TAlphabet>::VALUE;
    for (TSize i = 0; i < alphSize; ++i)
        std::cout << i << ',' << TAlphabet(i) << "  ";
    std::cout << std::endl;
}

int main()
{
	typedef String<AminoAcid> TAminoAcidString;
    TAminoAcidString str = "MQDRVKRPMNAFIVWSRDQRRKMALEN";
    Iterator<String<AminoAcid> >::Type it = begin(str);
	Iterator<String<AminoAcid> >::Type itEnd = end(str);

	std::cout << value(str);
	while(it != itEnd)
	{
		++it;
		if(value(it) == 'R')
		{
			std::cout<< 'A';
		} else {
		std::cout << *it;
		}
	}

	
	showAllLetterOfMyAlphabet(AminoAcid());
	
    showAllLetterOfMyAlphabet(Dna());
    showAllLetterOfMyAlphabet(Dna5());
    return 0;
}