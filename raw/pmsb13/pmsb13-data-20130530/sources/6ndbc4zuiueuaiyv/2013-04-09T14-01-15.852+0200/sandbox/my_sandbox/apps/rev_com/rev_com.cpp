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
    
    // Your code snippet
    resize(revComplGenome, sizeof(genome));
    std::cout << sizeof(revComplGenome);
    // And to check if your output is correct, 
    // use the given SeqAn function reverseComplement(),
    // which modifies the sequence in-place
    reverseComplement(genome);
    std::cout << genome << std::endl;
    return 0;
}

