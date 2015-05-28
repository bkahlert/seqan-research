#include <seqan/sequence.h>
#include <seqan/index.h>

using namespace seqan;

int main()
{
    String<Dna5> genome = "TTATTAAGCGTATAGCCCTATAAATATAA";
    StringSet<String<Dna5> > Patterns;
    resize(Patterns,3);
    Patterns[1]="TATAA";
    Patterns[2]="GCCC";
    Patterns[0]="CG";
    Index<String<Dna5>, IndexEsa<> > esaIndex(genome); 
    Finder<Index<String<Dna5>, IndexEsa<> > > esaFinder(esaIndex);
  
    for(int i=0;i<length(Patterns);i++){
      std::cout << Patterns[i] << '\n';
      while(find(esaFinder,Patterns[i] )){  
      std::cout << position(esaFinder) << '\t'; // -> 0
      }
      std::cout <<'\n'; // -> 0
      clear(esaFinder);
    }  
    return 0;
}