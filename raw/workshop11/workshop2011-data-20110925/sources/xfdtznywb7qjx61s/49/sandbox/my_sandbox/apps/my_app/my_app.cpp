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

template <typename String> class <CharString, DnaString, Peptide>
void countOneMers(String const& str)
{
	String<int> table;
	resize(table, 'z' - 'a' +1,0); 
	for(unsigned int i =0; i< length(str); i++)
	{
		table[str[i] - 'a'] += 1;
	}
	
	for (unsigned int i = 0; i < 'z' - 'a' + 1; ++i) {
        if (table[i] == 0)
            continue;
        std::cout << static_cast<char>('a' + i) << " " << table[i] << std::endl;
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