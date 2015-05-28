#include <iostream>
#include <seqan/sequence.h>
#include <seqan/file.h>

using namespace seqan;

int main()
{
    String<Dna5> nucleotides = "AGTCGTGNNANCT";
    String<Dna5> lesser;
    String<Dna5> greater;
    
    Dna5 splitNuc = 'G';
    for (unsigned i = 0; i < length(nucleotides); ++i){
	Lexical<> comp(nucleotides[i],splitNuc);
	if (isLess(comp))
	    appendValue(lesser, nucleotides[i]);
	if (isGreater(comp))
	    appendValue(greater, nucleotides[i]);
    }
    std::cout << "lesser: " << lesser << std::endl;
    std::cout << "greater: " << greater << std::endl;
    return 0;
}