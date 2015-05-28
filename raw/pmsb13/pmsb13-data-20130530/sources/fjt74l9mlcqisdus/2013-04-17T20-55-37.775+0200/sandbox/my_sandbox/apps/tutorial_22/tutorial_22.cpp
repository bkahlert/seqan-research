#include <seqan/seq_io.h>
#include <seqan/find_motif.h>
#include <seqan/index.h>

using namespace std;
using namespace seqan;


int main(){
	
	StringSet<CharString> reads;
	append(reads,"MHI");
	append(reads,"AMH");
	StringSet<CharString> protein_seq;
	append(protein_seq,"ILHIVAMHIIM");
	append(protein_seq,"IIMAMHIVAMI");
	
	
	Index<StringSet<CharString>> index_esa(protein_seq);
	Index<StringSet<CharString>> index(reads);
	Finder<Index<StringSet<CharString>>>finder_esa(index_esa);
	Pattern<Index<StringSet<CharString>>pattern(index);
	cout << "hit at ";
    while (find(finder_esa, pattern))
	cout << position(finder_esa) << " ";
	cout << endl;
	
	return 0;
}