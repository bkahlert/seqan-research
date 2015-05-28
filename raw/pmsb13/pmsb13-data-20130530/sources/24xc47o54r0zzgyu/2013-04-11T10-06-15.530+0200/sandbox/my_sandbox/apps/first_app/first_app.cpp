#include <seqan/sequence.h>
#include <seqan/index.h>

using namespace seqan;

int main()
{
    String<Dna> text = "TTATTAAGCGTATAGCCCTATAAATATAA";
	StringSet<DnaString> patterns;
	appendValue(patterns, "TATAA");
	appendValue(patterns, "AGCG");
	appendValue(patterns, "CCTA");
    
	Index<String<Dna> > index(text);
	Finder<Index<String<Dna> > > esafinder(index);

	typedef Iterator<StringSet<DnaString> >::Type SIterator;
	for(SIterator pattern = begin(patterns); pattern!=end(patterns);goNext(pattern))
	{
		std::cout << "Occurences of pattern " << *pattern << ":" << std::endl;
		while(find(esafinder, pattern))
			std::cout<<"match at position"<<position(esafinder)<<std::endl;
		clear(esafinder);
	}

    return 0;
}