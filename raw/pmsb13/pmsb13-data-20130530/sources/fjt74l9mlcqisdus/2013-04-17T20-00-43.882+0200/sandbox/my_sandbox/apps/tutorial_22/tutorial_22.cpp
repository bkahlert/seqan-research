#include <seqan/seq_io.h>
#include <seqan/find_motif.h>

using namespace std;
using namespace seqan;

int main(){
	String<String<Dna>> reads;
	append(reads,"ACGTAGACGATGAC");
	append(reads,"ACGCGATAGACGTC");
	

	return 0;
}