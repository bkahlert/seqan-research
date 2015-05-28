#include <seqan/sequence.h>
#include <seqan/index.h>

using namespace seqan;

int main()
{
    String<Dna5> genome = "TTATTAAGCGTATAGCCCTATAAATATAA";
    Index<String<Dna5>, IndexEsa<> > esaIndex(genome); 
    Finder<Index<String<Dna5>, IndexEsa<> > > esaFinder(esaIndex);


    while (find == true) {
        find(esaFinder, "TATAA"); // first occurrence of "ACGT"  
        std::cout << position(esaFinder) << std::endl; // -> 0
    }

    return 0;
}
