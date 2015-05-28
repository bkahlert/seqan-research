#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <iostream>

using namespace seqan;

// We want to define a function, which takes
// the alphabet type as an argument
template <typename TAlphabet>
void showAllLettersOfMyAlphabet(TAlphabet const &)
{
    unsigned abcSize = ValueSize<TAlphabet>::VALUE;
    for (int i=0; i < abcSize; i++){
    	std::cout << "Letter " << i << ": " << TAlphabet(i) << std::endl;
    }
}

int main()
{
    showAllLettersOfMyAlphabet(AminoAcid());
    showAllLettersOfMyAlphabet(Dna());
    showAllLettersOfMyAlphabet(Dna5());
    return 0;
}
