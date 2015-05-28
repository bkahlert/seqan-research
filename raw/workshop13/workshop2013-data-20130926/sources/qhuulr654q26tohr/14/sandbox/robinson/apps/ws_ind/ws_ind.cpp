#include <seqan/sequence.h>
#include <seqan/index.h>
#include <iostream>

using namespace std;
using namespace seqan;

int main()
{
    String<Dna> text = "TTATTAAGCGTATAGCCCTATAAATATAA";
    String<Dna> pattern = "TATAA";
    Index<String<Dna>, IndexEsa< > > esaIndex(text);
    Finder<Index<String<Dna>, IndexEsa<> > > esaFinder(esaIndex);

    while (find(esaFinder, pattern)) {
    	cout << "Found " << pattern << " at position " << position(esaFinder) << endl;
    }



    return 0;
}