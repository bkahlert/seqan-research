#include <seqan/seq_io.h>

seqan::FaiIndex faiIndex;
// ... index building here ...
int main(){
	int res = write(faiIndex, "path/to/file.fasta.fai");
	if (res != 0)
		std::cerr << "ERROR: Could not write the index to file!\n";
}
