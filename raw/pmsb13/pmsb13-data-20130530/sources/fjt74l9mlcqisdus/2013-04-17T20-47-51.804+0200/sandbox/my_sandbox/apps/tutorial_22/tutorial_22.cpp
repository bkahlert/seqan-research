#include <seqan/seq_io.h>
#include <seqan/find_motif.h>
#include <seqan/index.h>

using namespace std;
using namespace seqan;


int main(){
	
	String<CharString> reads;
	append(reads,"MHI");
	append(reads,"AMH");
	String<CharString> protein_seq;
	append(protein_seq,"ILHIVAMHIIM");
	append(protein_seq,"IIMAMHIVAMI");
	
	
	Index<String<CharString>> index_esa(protein_seq);
	Finder<Index<String<CharString>>>finder_esa(index_esa);
	Pattern<reads>pattern;
	cout << "hit at ";
    while (find(finder_esa, pattern))
	cout << position(finder_esa) << " ";
	cout << endl;
	
	return 0;
}