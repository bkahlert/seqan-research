#include <iostream>
#include <seqan/sequence.h>
#include <seqan/file.h>

using namespace seqan;

int main()
{
    String<Dna5> nucleotides = "AGTCGTGNNANCT";
    String<Dna5> lesser;
	String<Dna5> greater;

    // Append all elements of nucleotides, apart of Gs, 
    // to the list `selected` 
	for (unsigned i = 0; i < length(nucleotides); ++i){
        if (nucleotides[i]<(Dna5)'G') appendValue(lesser, nucleotides[i]);
		if (nucleotides[i]>(Dna5)'G') appendValue(greater, nucleotides[i]);
    }
    std::cout << "lesser nucleotides: " << lesser << std::endl;
    std::cout << "greater nucleotides: " << greater << std::endl;
    system("Pause");
	return 0;
}