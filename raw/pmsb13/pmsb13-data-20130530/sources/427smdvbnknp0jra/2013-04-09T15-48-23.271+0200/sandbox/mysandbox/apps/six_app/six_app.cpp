#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <iostream>

using namespace seqan;

// We want to define a function, which takes 
// the alphabet type as an argument
template <typename TAlphabet>
void showAllLettersOfMyAlphabet(TAlphabet const & input)
{
	
	unsigned alphSize = ValueSize<TAlphabet>::VALUE;
	std::cout << "Alphabet size : \n" << alphSize << std::endl;   
    for (unsigned i = 0; i < alphSize; ++i)
        std::cout << i << ': ' << TAlphabet(i) << "\n";
    std::cout << std::endl;
}

int main()
{
    showAllLettersOfMyAlphabet(AminoAcid());
    showAllLettersOfMyAlphabet(Dna());
    showAllLettersOfMyAlphabet(Dna5());
    return 0;
}
