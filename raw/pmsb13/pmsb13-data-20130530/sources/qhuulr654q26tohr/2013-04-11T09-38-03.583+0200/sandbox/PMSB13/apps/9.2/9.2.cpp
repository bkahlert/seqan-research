#include <seqan/sequence.h>
#include <seqan/index.h>

using namespace seqan;
using namespace std;

int main()
{
    String<Dna> haystack = "TTATTAAGCGTATAGCCCTATAAATATAA";
    String<Dna> needle = "TATAA";
    typedef Index<String<Dna>, FMIndex<> > TIndex;
    typedef Finder<TIndex> TFinder;
    TIndex index(haystack);
    TFinder finder(index);
    
    cout << "Found occurences at:\t";
    while (find(finder,needle))
	cout << position(finder) << " ";
    cout << endl;
    
    return 0;
}