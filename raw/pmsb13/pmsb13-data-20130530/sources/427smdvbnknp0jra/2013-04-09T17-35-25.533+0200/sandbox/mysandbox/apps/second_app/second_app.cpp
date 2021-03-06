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
	std::cout << genome << std::endl;
    Iterator<DnaString, Rooted >::Type it = begin(genome);
	ModifiedString< DnaString, ModReverse > myModifier(genome);
	std::cout << *it << std::endl;
	for (; !atEnd(it);it++)
		*it=	getRevCompl(*it);
	

	std::cout << myModifier << std::endl;
	
    //reverseComplement(genome);
    
    return 0;
}