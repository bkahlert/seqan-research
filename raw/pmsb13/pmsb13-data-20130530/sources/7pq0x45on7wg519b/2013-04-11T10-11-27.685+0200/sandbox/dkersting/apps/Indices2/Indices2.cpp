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
	StringSet<String<Dna5> > pattern;
	appendValue(pattern,  "TATA");
	appendValue(pattern,  "GCGT");
	appendValue(pattern,  "TAGC");

	Index<StringSet<String<Dna5> >, FMIndex> fmIndex(pattern); 
    	Finder<Index<String<Dna5>, IndexEsa<> > > esaFinder(esaIndex);

    typedef Iterator<StringSet>::Type TIterator;
	for (TIterator it = begin(mySet); it != end(mySet); goNext(it))
	{
		std::cout << "Find pattern "<<getValue(it)<<" at Position:"<<std::endl;
 		//while(find(esaFinder, ))
    		//{
        	//std::cout << position(esaFinder) << std::endl;
    		//}
	}

}
