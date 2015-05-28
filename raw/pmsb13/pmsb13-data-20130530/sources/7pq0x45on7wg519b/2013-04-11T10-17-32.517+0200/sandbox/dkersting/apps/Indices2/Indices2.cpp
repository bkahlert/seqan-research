// ==========================================================================
//                                  Indices2
// ==========================================================================
#include <seqan/sequence.h>
#include <seqan/index.h>

using namespace seqan;
int main()
{
	String<Dna5> genome = "TTATTAAGCGTATAGCCCTATAAATATAA";
	Index<String<Dna5>, IndexEsa<> > esaIndex(genome);
	Finder<Index<String<Dna5>, IndexEsa<> > > esaFinder(esaIndex);
	StringSet<String<Dna5> > pattern;
	appendValue(pattern,  "TATA");
	appendValue(pattern,  "GCGT");
	appendValue(pattern,  "TAGC");

    typedef Iterator<StringSet<String<Dna5> > >::Type TIterator;
	for (TIterator it = begin(pattern); it != end(pattern); goNext(it))
	{
		std::cout << "Find pattern "<<getValue(it)<<" at Position:"<<std::endl;
 		while(find(esaFinder, ))
    		{
        		std::cout << position(esaFinder) << std::endl;
    		}
		clear(esaFinder);
	}

}
