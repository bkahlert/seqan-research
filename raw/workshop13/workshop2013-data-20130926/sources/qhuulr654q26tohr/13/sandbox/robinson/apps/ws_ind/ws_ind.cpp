#include <seqan/sequence.h>
#include <seqan/index.h>
#include <iostream>

using namespace std;
using namespace seqan;

int main()
{
    String<char> text = "TTATTAAGCGTATAGCCCTATAAATATAA";
    String<char> pattern = "TATAA";
    Index<String<char>, IndexEsa< > > esaIndex(text);
    Finder<Index<String<char>, IndexEsa<> > > esaFinder(esaIndex);

    while (find(esaFinder, pattern)) {
    	cout << "Found " << pattern << " at position " << position(esaFinder) << endl;
    }



    return 0;
}