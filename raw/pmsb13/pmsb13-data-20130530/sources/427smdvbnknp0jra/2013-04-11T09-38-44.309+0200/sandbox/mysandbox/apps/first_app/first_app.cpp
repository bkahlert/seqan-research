#include <seqan/sequence.h>
#include <seqan/index.h>

using namespace seqan;

int main()
{
    String<Dna> gen = "TTATTAAGCGTATAGCCCTATAAATATAA";
	String<Dna> pattern = "TATAA";

    Index<String<Dna>, FMIndex<> > index(gen);
	Finder<Index<String<Dna>, FMIndex<> >> FMFinder(index);

	while(find(FMFinder,pattern))
	{
		std::cout << position(FMFinder) << std::endl;
	
	}

    return 0;
}