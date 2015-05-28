#include <seqan/sequence.h>
#include <seqan/index.h>

using namespace seqan;

int main()
{
	String<char> text = "TTATTAAGCGTATAGCCCTATAAATATAA";    
	String<char> pattern = "TATAA";

	Index<String<char> > esaIndex(text);
	Findex<Index<String<char> > > esaFinder(esaIndex);
	
	while(find(esaFinder, pattern))
		std::cout << position(esaFinder) <<std::endl;

    return 0;
}
