#include <seqan/sequence.h>
#include <seqan/index.h>

using namespace seqan;

int main()
{
    String<Dna5> genome = "ACGTACGTACGTN";
    Index<String<Dna5>, IndexEsa > esaIndex(genome); 
    Finder<Index<String<Dna5>, IndexEsa> > esaFinder(esaIndex);

    find(esaFinder, "ACGT"); // first occurrence of "ACGT"  
    position(esaFinder); // -> 0
    find(esaFinder, "ACGT"); // second occurrence of "ACGT"
    position(esaFinder); // -> 4
    find(esaFinder, "ACGT"); // third occurrence of "ACGT"
    position(esaFinder); // -> 8
}
