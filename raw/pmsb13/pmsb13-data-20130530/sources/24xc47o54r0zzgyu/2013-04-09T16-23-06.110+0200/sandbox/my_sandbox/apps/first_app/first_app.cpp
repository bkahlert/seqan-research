#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <iostream>

using namespace seqan;
using namespace std;

// We want to define a function, which takes 
// the alphabet type as an argument
template <typename TAlphabet>
void showAllLettersOfMyAlphabet(TAlphabet const &)
{
    unsigned alphSize = ValueSize<TAlphabet>::VALUE;
	for (unsigned i = 0; i<alphSize;++i)
	{
		cout << TAlphabet(i)<< endl;
	}
}

int main()
{
    showAllLettersOfMyAlphabet(AminoAcid());
    showAllLettersOfMyAlphabet(Dna());
    showAllLettersOfMyAlphabet(Dna5());
    return 0;
}