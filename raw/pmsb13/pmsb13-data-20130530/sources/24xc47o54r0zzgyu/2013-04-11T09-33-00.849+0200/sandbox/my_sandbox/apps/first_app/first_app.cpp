#include <seqan/sequence.h>
#include <seqan/index.h>

using namespace seqan;

int main()
{
    String<Dna> text = "TTATTAAGCGTATAGCCCTATAAATATAA";
	String<Dna> pattern = "TATAA";
    
	Index<String<Dna> > index(text);
	Finder<Index<String<Dna> > > esafinder(index);

	while(find(esafinder, pattern))
		std::cout<<position(esafinder)<<std::endl;

    return 0;
}