#include <seqan/seq_io.h>
#include <seqan/stream.h>

using namespace std;

int main(int argc,char const ** argv){

	if (argc>=2){
		seqan::FaiIndex faiIndex;
		int res = read(faiIndex, "D:\\fasta1.fa");
		if (res != 0) cerr << "ERROR: Could not read FAI index path/to/file.fasta.fai\n";
		unsigned idx = 0;
		if (!getIdByName(faiIndex, argv[2], idx)) cerr << "ERROR: FAI index has no entry for chr1.\n";
		
		unsigned seqLength = sequenceLength(faiIndex, idx);

		// Load first 1000 characters of chr1.
		seqan::CharString seqinfix;
		int begin,end;
		seqan::lexicalCast2(begin, argv[3]);
		seqan::lexicalCast2(end, argv[4]);
		if (readRegion(seqinfix, faiIndex, idx, begin,end) != 0) cerr << "ERROR: Could not load chr1.\n";

		// Load all of chr1.
		seqan::CharString seqChr1;
		if (readSequence(seqChr1, faiIndex, idx) != 0) cerr << "ERROR: Could not load chr1.\n";
		
		
		cout << seqinfix<<seqan::length(seqinfix)<<endl;
	}
	else return 1;
	return 0;	
}



