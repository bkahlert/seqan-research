#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <iostream>
#include <seqan/file.h>
#include <seqan/modifier.h>

using namespace seqan;

Dna getRevCompl(Dna const & nucleotide)
{
    if(nucleotide == (Dna)'A')
        return (Dna)'T';
    if(nucleotide == (Dna)'T')
        return (Dna)'A';
    if(nucleotide == (Dna)'C')
        return (Dna)'G';
    return (Dna)'C';
}


int main()
{
    DnaString genome = "TATATACGCGCGAGTCGT";
    DnaString revComplGenome;
    resize(revComplGenome,length(genome));
	for (unsigned i=0;i < length(genome);i++)
		revComplGenome[length(genome)-i-1]=	getRevCompl(genome[i]);
	
	std::cout << genome << "\n" ;
	for (unsigned i = length(revComplGenome); i > 0; i--)
        std::cout << revComplGenome[i];
	std::cout << "\n\n" ;
	

    reverseComplement(genome);
    std::cout << genome << std::endl;
    return 0;
}