#include <seqan/sequence.h>
#include <seqan/index.h>
#include <iostream>

using namespace std;
using namespace seqan;

int main()
{
    String<char> text = "This is the first example";
    Index<String<char>, IndexEsa<> > index(text);

    String<Dna> genome= "TTATTAAGCGTATAGCCCTATAAATATAA";
    Index<String<Dna>, IndexEsa<> > esaIndex(genome); 
    Finder<Index<String<Dna>, IndexEsa<> > > esaFinder(esaIndex);
    
   
   while find(esaFinder, "TATAA"){
     cout << position(esaFinder) << endl;
   }
    return 0;
}