#include <seqan/sequence.h>
#include <seqan/index.h>
#include <fstream>
#include <string>

using namespace seqan;

int main()
{
    String<Dna> gen = "TTATTAAGCGTATAGCCCTATAAATATAA";
	String<Dna> pattern = "TATAA";

	//"TATAA\nTAA\n>CCT\nAGC\nTT";
	
	std::cout << "FMIndex" << std::endl;

    Index<String<Dna>, FMIndex<> > index(gen);
	Finder<Index<String<Dna>, FMIndex<> >> FMFinder(index);

	while(find(FMFinder,pattern))
	{
		std::cout << position(FMFinder) << std::endl;
	
	}

	std::cout << "\n" << "EsaIndex" << std::endl;
	pattern="CCT";
	Index<String<Dna>, IndexEsa<> > esaIndex(gen);
	Finder<Index<String<Dna>, IndexEsa<> >> esaFinder(esaIndex);

	while(find(esaFinder,pattern))
	{
		std::cout << position(esaFinder) << std::endl;
	
	}

    return 0;
}