#include <seqan/sequence.h>
#include <seqan/index.h>

using namespace seqan;

int main()
{
    String<Dna> dna = "TTATTAAGCGTATAGCCCTATAAATATAA";
	String<Dna> pattern = "TATAA";
    Index<String<Dna>, IndexEsa<> > index(dna);
	Finder<Index<String<Dna>>, IndexEsa<>> esaFinder(index);
	while(find(esaFinder,"TATAA")){
		std::cout << position(esaFinder) << std::endl;}


	/*
	StringSet<String<Dna> > dna2;
	appendValue(dna2,  "ACGT");
	appendValue(dna2,  "CTAG");
	appendValue(dna2,  "GAAA");

	Index<StringSet<String<Dna> >, IndexEsa<>> fmIndex(dna2); 
	*/


    return 0;
}