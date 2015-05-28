#include <iostream>
#include <seqan/sequence.h>
#include <seqan/file.h>

using namespace seqan;

int main()
{
    String<Dna5> nucleotides = "AGTCGTGNNANCT";
	String<Dna5> g="G";
	String<Dna5> lesser;
	String<Dna5> greater;
    // Append all elements of nucleotides, apart of Gs, 
    // to the list `selected` 
    for (unsigned i = 0; i < length(nucleotides); ++i){
		if(nucleotides[i] < g)
			appendValue(lesser, nucleotides[i]);
		else
			if(nucleotides[i] != g)
				appendValue(greater, nucleotides[i]);
    }
    std::cout << "All nucleotides: " << nucleotides << "\n" << "Unselected nucleotides: " << g << "\n" << "Lesser nucleotides: "<< lesser << "\n" << "Greater nucleotides: " << greater << std::endl;
    return 0;
}