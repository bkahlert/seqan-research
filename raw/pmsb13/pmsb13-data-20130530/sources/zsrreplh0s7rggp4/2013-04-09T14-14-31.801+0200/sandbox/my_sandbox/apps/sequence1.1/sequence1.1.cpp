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




DnaString revCompl(DnaString const & genome)
  {
   DnaString RevCompl;
   seqan::resize(RevCompl,length(genome));
   for (int i=1,length(genome),i++){
    RevCompl[length(genome)-1]=getRevCompl(genome[i-1]);
  }

  



int main()
{
  
  revCompl(genome);
  DnaString genome = "TATATACGCGCGAGTCGT";
    
    
    
    
    DnaString revComplGenome;
    
    // Your code snippet

    // And to check if your output is correct, 
    // use the given SeqAn function reverseComplement(),
    // which modifies the sequence in-place
   // reverseComplement(genome);
    std::cout << genome==RevCompl << std::endl;
    return 0;
}