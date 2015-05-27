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
	String<char> str = "MQDRVKRPMNAFIVWSRDQRRKMALEN";
    Iterator<String<char> >::Type it = begin(str);
	Iterator<String<char> >::Type itEnd = end(str);
while(it != itEnd)
{
	std::cout << *it;
	
	if(*it == "M")
	{
	std::cout << "\n";	
	}
	++it;
}
	
	showAllLetterOfMyAlphabet(AminoAcid());
	
    showAllLetterOfMyAlphabet(Dna());
    showAllLetterOfMyAlphabet(Dna5());
    return 0;
}