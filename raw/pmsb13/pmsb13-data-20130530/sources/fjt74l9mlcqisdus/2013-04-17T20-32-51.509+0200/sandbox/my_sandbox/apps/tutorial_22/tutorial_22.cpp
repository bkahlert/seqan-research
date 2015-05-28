#include <seqan/seq_io.h>
#include <seqan/find_motif.h>
#include <seqan/index.h>

using namespace std;
using namespace seqan;


int main(){
	String<String<char>> reads;
	append(reads,"MHI");
	append(reads,"AMH");
	String<String<char>> protein_seq;
	append(protein_seq,"ILHIVAMHIIM");
	append(protein_seq,"IIMAMHIVAMI");
	
	
	Index<String<char>> index_esa(protein_seq);
	Finder<Index<String<char>>>finder_esa(index_esa);

	cout << "hit at ";
    while (find(finder_esa, reads))
	cout << position(finder_esa) << " ";
	cout << endl;
	return 0;
}