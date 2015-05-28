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
    
    for (unsigned i = 0; i < length(patternSet); i++) {
        std::cout << "Positions for " << patternSet[i] << ":\t";
        while (find(esaFinder, patternSet[i]) == 1) {
            std::cout << position(esaFinder) << "\t";
        }
        clear(esaFinder);
        std::endl;
    }

    return 0;
}
