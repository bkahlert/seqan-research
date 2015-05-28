#include <seqan/sequence.h>
#include <seqan/index.h>

using namespace seqan;
using namespace std;

int main()
{
    String<Dna> genome = "TTATTAAGCGTATAGCCCTATAAATATAA";
	String<Dna> pattern = "TATAA";
	Index<String<Dna>> index(genome);
	Finder<Index<String<Dna>>> defaultfinder(index);
	
	while(find(defaultfinder,pattern)){
		cout << position(defaultfinder)<<endl;
	}
	
    return 0;
}