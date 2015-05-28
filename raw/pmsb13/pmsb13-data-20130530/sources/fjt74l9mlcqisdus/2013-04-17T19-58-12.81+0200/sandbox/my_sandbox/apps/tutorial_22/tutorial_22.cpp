#include <seqan/seq_io.h>
#include <seqan/find_motif.h>

using namespace std;
using namespace seqan;

int main(){
	String<String<Dna>> reads;
	addValue(reads,"ACGTAGACGATGAC");
	addValue(reads,"ACGCGATAGACGTC");
	

	return 0;
}