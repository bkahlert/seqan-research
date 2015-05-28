#include <seqan/sequence.h>
#include <seqan/index.h>

using namespace seqan;

int main()
{
	
    String<Dna5> dna = "TTATTAAGCGTATAGCCCTATAAATATAA";
	String<Dna5> pattern = "TATAA";
    Index<String<Dna5>, IndexEsa<> > index(dna);
	Finder<Index<String<Dna5>, IndexEsa<> >> esaFinder(index);
	while(find(esaFinder,"TATAA")){
		std::cout << position(esaFinder) << std::endl;}


	/*
	StringSet<String<Dna> > dna2;
	appendValue(dna2,  "ACGT");
	appendValue(dna2,  "CTAG");
	appendValue(dna2,  "GAAA");

	Index<StringSet<String<Dna> >, IndexEsa<>> fmIndex(dna2); 

	/*
    String<Dna5> genome = "TTATTAAGCGTATAGCCCTATAAATATAA";
    Index<String<Dna5>, IndexEsa<> > esaIndex(genome);
    Finder<Index<String<Dna5>, IndexEsa<> > > esaFinder(esaIndex);
    while(find(esaFinder, "TATAA"))
    {
        std::cout << position(esaFinder) << std::endl;
    }*/
    return 0;
}


    return 0;
}