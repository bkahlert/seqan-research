#include <iostream>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <string>
#include <seqan/stream.h>

int main(int argc,char *argv[]){
	if(argc <5){
		std::cerr << "more input!";
		return 1;
	}

	// Try to load index and create on the fly if necessary.
    seqan::FaiIndex faiIndex;
    if (seqan::read(faiIndex, argv[1]) != 0)
    {
        std::cerr << "ERROR: Index could not be loaded.\n";
		return 1;
    }
    // Translate sequence name to index.
    unsigned idx = 0;
    if (!getIdByName(faiIndex, argv[2], idx))
    {
        std::cerr << "ERROR: Index does not know about sequence " << argv[2] << "\n";
        return 1;
    }
    // Convert positions into integers.
    unsigned beginPos = int(argv[3]);
	unsigned endPos = int(argv[4]);
    
    // Make sure begin and end pos are on the sequence and begin <= end.
    if (beginPos > sequenceLength(faiIndex, idx))
        beginPos = sequenceLength(faiIndex, idx);
    if (endPos > sequenceLength(faiIndex, idx))
        endPos = sequenceLength(faiIndex, idx);
    if (beginPos > endPos)
        endPos = beginPos;
    // Finally, get infix of sequence.
    seqan::Dna5String sequenceInfix;
    if (readRegion(sequenceInfix, faiIndex, idx, beginPos, endPos) != 0)
    {
        std::cerr << "ERROR: Could not load infix.\n";
        return 1;
    }
    std::cout << sequenceInfix << "\n"
	


	return 0;
}