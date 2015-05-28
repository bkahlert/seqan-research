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

void my_reverseComplement(DnaString & referenz, DnaString & revcomp){
	int len = seqan::length(referenz)-1;
	for (int i= 0;i<=len;i++){
		revcomp[len-i]=getRevCompl(referenz[i]);
	}
}

int main()
{
    DnaString genome = "TATATACGCGCGAGTCGT";
    DnaString revComplGenome;
	resize(revComplGenome,seqan::length(genome));
    my_reverseComplement(genome,revComplGenome);
    //reverseComplement(genome);
    std::cout << genome << std::endl;
	std::cout << revComplGenome << std::endl;
    system("Pause");
	return 0;
}