#include <seqan/sequence.h>
#include <seqan/index.h>

using namespace seqan;
using namespace std;

int main()
{
    String<Dna> haystack = "TTATTAAGCGTATAGCCCTATAAATATAA";
    typedef StringSet<String<Dna> > TNeedles;
    TNeedles needles;
    appendValue(needles,"TATAA");
    appendValue(needles,"TAA");
    appendValue(needles,"GC");
    typedef Index<String<Dna>, FMIndex<> > TIndex;
    typedef Finder<TIndex> TFinder;
    TIndex index(haystack);
    TFinder finder(index);
    
    typedef Iterator<TNeedles, Rooted> TIteratorNeedle;
    TIteratorNeedle it = begin(needles);
    
    for (;atEnd(it);goNext(it)){
	clear(finder);
	cout << "Found occurences of "<<*it<<" at:\t";
	while (find(finder,*it))
	    cout << position(finder) << " ";
	cout << endl;
    }
    
    
    return 0;
}