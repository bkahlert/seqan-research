#include <iostream>
#include <seqan/sequence.h>
#include <seqan/file.h>

using namespace seqan;

int main()
{
    String<Dna5> nucleotides = "AGTCGTGNNANCT";
    String<Dna5> lesser;
	String<Dna5> greater; 
    for (unsigned i = 0; i < length(nucleotides); ++i)
	{
		if(nucleotids[i] < 'G')
        	appendValue(lesser, nucleotides[i]);
		else
			appendValue(greater, nucleotides[i]);
    }
	std::cout << "Lexicographically smaller: " << lesser << std::endl;  
	std::cout << "Lexicographically greater: " << greater << std::endl;
    return 0;
}
