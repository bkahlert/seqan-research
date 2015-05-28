#include <iostream>
#include <seqan/sequence.h>
#include <seqan/file.h>

using namespace seqan;

int main()
{
    String<Dna5> nucleotides = "AGTCGTGNNANCT";
    String<Dna5> selected;
    String<Dna5> smaller;
    String<Dna5> greater;
    // Append all elements of nucleotides, apart of Gs,
    // to the list `selected`
    for (unsigned i = 0; i < length(nucleotides); ++i){
        if (nucleotides[i] < 'G')
            appendValue(smaller, nucleotides[i]);
        else if (nucleotides[i] > 'G')
            appendValue(greater, nucleotides[i]);
    }
    std::cout << "smaller nucleotides: " << smaller << std::endl;
    std::cout << "greater nucleotides: " << greater << std::endl;
    return 0;
}
