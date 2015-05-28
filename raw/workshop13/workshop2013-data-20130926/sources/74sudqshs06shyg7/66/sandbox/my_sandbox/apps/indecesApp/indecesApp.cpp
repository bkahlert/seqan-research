#include <seqan/string.h>
#include <seqan/index.h>

using namespace seqan;

int main()
{
    String<Dna5> genome = "TTATTAAGCGTATAGCCCTATAAATATAA";
    Index<String<Dna5>, FirstIndex> DnaIndex(genome); 
    Finder< Index<String<Dna5>, IndexEsa> > DnaFinder(DnaIndex);

    find(DnaFinder, "TATA");
    int pos=position(DnaFinder);
    
    std::cout << pos << std::endl;
    
}