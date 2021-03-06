#include <seqan/sequence.h>
#include <seqan/index.h>

using namespace seqan;

int main()
{
    String<Dna5> genome = "TTATTAAGCGTATAGCCCTATAAATATAA";
    Index< String<Dna5>, IndexEsa<> > DnaIndex(genome); 
    Finder< Index<String<Dna5>, IndexEsa<> > > DnaFinder(DnaIndex);

    while(find(DnaFinder, "TATA")){
    
    	int pos=position(DnaFinder);
    
   	std::cout << pos << std::endl;
   }
    
}