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
	typedef ModifiedString<String<Dna>, ModComplementDna>   TMyComplement;
	typedef ModifiedString<TMyComplement, ModReverse>       TMyReverseComplement;    
    ::std::cout << myReverseComplement << ::std::endl;
    
    reverseComplement(genome);
    std::cout << genome << std::endl;
    return 0;
}