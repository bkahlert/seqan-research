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
DnaString revComplement(DnaString const & gene)
{
	DnaString revComplGene;
	for(unsigned i = 0; i < length(gene); ++i)
	{ 	
		revComplGene[i] = getRevCompl(gene[i]);
	}	
	return revComplGene;
}
int main()
{
    DnaString genome = "TATATACGCGCGAGTCGT";
    DnaString revComplGenome;
    rezise(revComplGenome,length(genome));

	revComplGenome = revComplement(genome);
	std::cout << revComplGenome << std::endl;
    reverseComplement(genome);
    std::cout << genome << std::endl;
    return 0;
}
