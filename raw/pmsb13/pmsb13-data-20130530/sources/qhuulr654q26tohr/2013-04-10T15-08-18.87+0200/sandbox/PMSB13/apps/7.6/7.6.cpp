#include <iostream>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

using namespace seqan;
using namespace std;

int main(int argc, char** argv )
{
    if (argc!=2){
	std::cerr << "USAGE: 7.5 FILENAME\n";
	return 1;
    }
    
    
    seqan::SequenceStream seqStream(argv[1], seqan::SequenceStream::WRITE);
    
    if (!isGood(seqStream))
    {
	std::cerr << "ERROR: Could not open the file.\n";
	return 1;
    }
    
    
    seqan::StringSet<CharString> ids;
    seqan::StringSet<Dna5String> seqs;
    seqan::StringSet<CharString> quals;
    
    resize(ids,2);
    resize(seqs,2);
    
    appendValue(ids,"seq1");
    appendValue(seqs,"CGAT");
    appendValue(ids,"seq2");
    appendValue(seqs,"ACGT");
    
    if (writeAll(seqStream, ids, seqs) != 0)
	{
	    std::cerr << "ERROR: Could not write to " << argv[1] << "!\n";
	    return 1;
	}
	return 0;
}