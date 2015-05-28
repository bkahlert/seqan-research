#include <seqan/sequence.h>
#include <seqan/index.h>

using namespace seqan;

int main()
{
    String<Dna5> genome = "TTATTAAGCGTATAGCCCTATAAATATAA";
    Index<String<Dna5>, IndexEsa<> > esaIndex(genome); 
    Finder<Index<String<Dna5>, IndexEsa<> > > esaFinder(esaIndex);

    StringSet< String<char> > patternSet;
    resize(patternSet,3);
    patternSet[0] = "TATAA";
    patternSet[1] = "GATA";
    patternSet[2] = "AGCC";
    
    std::cout << length(patternSet) << std::endl;

/*    for (unsigned i = 0; i < sizeof(patternSet); i++) {
        while (find(esaFinder, patternSet[i]) == 1) {
            std::cout << position(esaFinder) << std::endl;
        }
    }
*/

    return 0;
}
