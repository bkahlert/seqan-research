#include <seqan/seq_io.h>


using namespace std;
using namespace seqan;

int main(){
	String<String<Dna>> reads;
	addValue(reads,"ACGTAGACGATGAC");
	addValue(reads,"ACGCGATAGACGTC");
	

	return 0;
}