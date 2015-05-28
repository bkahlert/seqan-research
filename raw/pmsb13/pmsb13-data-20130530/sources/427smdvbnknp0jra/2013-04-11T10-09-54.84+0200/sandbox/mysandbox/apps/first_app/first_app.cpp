#include <seqan/sequence.h>
#include <seqan/index.h>
#include <fstream>
#include <string>

using namespace seqan;

int main()
{
    String<Dna> gen = "TTATTAAGCGTATAGCCCTATAAATATAA";
	String<Dna> pattern = "TATAA";

	String<String<Dna>> patterns;
	append(patterns,"TATAA");
	append(patterns,"TAA");
	append(patterns,"CCT");
	
	
	std::cout << "FMIndex" << std::endl;

    Index<String<Dna>, FMIndex<> > index(gen);
	Finder<Index<String<Dna>, FMIndex<> >> FMFinder(index);
	for(unsigned i=0;i<length(patterns);i++)
	{

		while(find(FMFinder,patterns[i]))
		{
			std::cout << position(FMFinder) << std::endl;
		
		}
		clear(FMFinder);
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