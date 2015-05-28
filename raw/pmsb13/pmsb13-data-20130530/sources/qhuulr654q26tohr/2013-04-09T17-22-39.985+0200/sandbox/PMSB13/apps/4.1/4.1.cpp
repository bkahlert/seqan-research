#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <iostream>

using namespace seqan;
using namespace std;

// We want to define a function, which takes 
// the alphabet type as an argument
template <typename TAlphabet>
void showAllLettersOfMyAlphabet(TAlphabet const & a)
{
    /*for (int i=0; i<ValueSize<TAlphabet>::VALUE;++i){
	std::cout << ordValue(TAlphabet(i)) << endl;
    }*/
}

int main()
{
    showAllLettersOfMyAlphabet(AminoAcid());
    showAllLettersOfMyAlphabet(Dna());
    showAllLettersOfMyAlphabet(Dna5());
    std::cout << Dna(6) << std::endl;
    return 0;
}