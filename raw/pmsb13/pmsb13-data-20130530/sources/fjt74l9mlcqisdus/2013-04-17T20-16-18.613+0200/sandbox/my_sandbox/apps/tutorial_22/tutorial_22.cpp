#include <seqan/seq_io.h>
#include <seqan/find_motif.h>

using namespace std;
using namespace seqan;

int main(){
	String<String<Dna>> reads;
	append(reads,"ACGTAGACGATGAC");
	append(reads,"ACGCGATAGACGTC");
	String<String<AminoAcid>> protein_seq;
	append(protein_seq,"ILHIVAMHIIM");
	append(protein_seq,"IIMAMHIVAMI");
	cout<<reads<<endl;
	cout<<protein_seq<<endl;

	return 0;
}