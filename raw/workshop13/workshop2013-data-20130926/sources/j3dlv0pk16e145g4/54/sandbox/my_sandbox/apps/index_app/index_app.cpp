#include <seqan/sequence.h>
#include <seqan/index.h>

using namespace seqan;

int main()
{
    String<Dna> text = "TTATTAAGCGTATAGCCCTATAAATATAA";
    String<Dna> pattern = "TATAA";
    Index<String<Dna> > index(text);
    Finder<Index<String<Dna> > > finder(index);
    while(find(finder, pattern))
    {
        std::cout << "occurence at position: " << position(finder) << std::endl;
    }
    return 0;
}