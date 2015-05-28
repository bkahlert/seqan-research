#include <seqan/seq_io.h>
#include <seqan/find_motif.h>
#include <seqan/index.h>

using namespace std;
using namespace seqan;


int main(){
	
	CharString reads= "HII";
	StringSet<CharString> protein_seq;
	resize(protein_seq,2);
	protein_seq[0]="ILHIVAMHIIM";
	protein_seq[1]="IIHIIHIVAMI";
	
	
	Index<StringSet<CharString>> index_esa(protein_seq);
	Finder<Index<StringSet<CharString>>>finder_esa(index_esa);
	Pattern<CharString>pattern(reads);
	cout << "hit at ";
    while (find(finder_esa, pattern))
	cout << position(finder_esa) << " ";
	cout << endl;
	
	return 0;
}