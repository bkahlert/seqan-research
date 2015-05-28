#include <seqan/seq_io.h>
#include <seqan/find_motif.h>
#include <seqan/index.h>

using namespace std;
using namespace seqan;


int main(){
	
	StringSet<CharString> reads;
	resize(reads,2);
	reads[0]="MHI";
	reads[1]="AMH";
	StringSet<CharString> protein_seq;
	resize(protein_seq,2);
	protein_seq[0]="ILHIVAMHIIM";
	protein_seq[1]="IIMAMHIVAMI";
	
	
	Index<StringSet<CharString>> index_esa(protein_seq);
	//Index<StringSet<CharString>> index(reads);
	Finder<Index<StringSet<CharString>>>finder_esa(index_esa);
	Pattern<StringSet<CharString>>pattern(reads);
	cout << "hit at ";
    while (find(finder_esa, pattern))
	cout << position(finder_esa) << " ";
	cout << endl;
	
	return 0;
}