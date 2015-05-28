#include <seqan/sequence.h>
#include <seqan/index.h>

using namespace seqan;
using namespace std;

int main()
{
    String<Dna> genome = "TTATTAAGCGTATAGCCCTATAAATATAA";
	StringSet<String<Dna>> pattern;
	appendValue(pattern,"TATAA");
	appendValue(pattern,"AAGCG");
	appendValue(pattern,"ATAGC");
	Index<String<Dna>> index(genome);
	Finder<Index<String<Dna>>> defaultfinder(index);
	for (int i=0;i<length(pattern);++i){
		cout << pattern[i]<<": ";
		while(find(defaultfinder,pattern[i])){
			cout << position(defaultfinder);
		}
		clear(defaultfinder);
		cout<<endl;
	}
    return 0;
}