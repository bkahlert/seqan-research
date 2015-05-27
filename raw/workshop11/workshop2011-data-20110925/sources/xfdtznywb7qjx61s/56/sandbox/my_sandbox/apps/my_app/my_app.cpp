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

//template <typename String> class <CharString, DnaString, Peptide>
template <typename TString>
void countOneMers(TString const sequence)
{
	typedef typename Size<TString>::Type TSize;
	typedef typename Value<TString>::Type TAlphabet;
	
	typedef String<TSize> TCounterString;
	TSize alphSize = ValueSize<TAlphabet>::VALUE;
	TCounterString table;
	resize(table, alphSize,0); 
	
	typedef typename Iterator<TString>::Type TIter;
    TIter itEnd = end(sequence);
    for (TIter it = begin(sequence); it != itEnd; goNext(it))
        value(table, ordValue(value(it))) += 1;
	
	
	typedef typename Iterator<TCounterString>::Type TTableIter;	
	TTableIter countIt = begin(table);
    TTableIter countItEnd = end(table);
	for(TSize pos=0; countIt != countItEnd; ++countIt, ++pos)
	{
		if(value(countIt > 0)
		{
			std::cout <<TAlphabet(pos) << ':' << value(countIt) << std::endl;
		}
	}

}
	

int main()
{
	countOneMers("mississipi");
	
	typedef String<AminoAcid> TAminoAcidString;
    TAminoAcidString str = "MQDRVKRPMNAFIVWSRDQRRKMALEN";
    Iterator<String<AminoAcid> >::Type it = begin(str);
	Iterator<String<AminoAcid> >::Type itEnd = end(str);

//	std::cout << str;
	while(it != itEnd)
	{
		if(value(it) == 'R')
		{
			std::cout<< 'A';
		} else {
		std::cout << *it;
		}
		++it;

	}

	
	showAllLetterOfMyAlphabet(AminoAcid());
    showAllLetterOfMyAlphabet(Dna());
    showAllLetterOfMyAlphabet(Dna5());
    return 0;
}