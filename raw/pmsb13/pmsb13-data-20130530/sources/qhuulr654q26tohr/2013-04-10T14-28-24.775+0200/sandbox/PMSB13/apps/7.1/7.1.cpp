#include <iostream>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

int main(int argc, char** argv )
{
    if (argc>2)
	return 1;
    
    seqan::CharString id;
    seqan::Dna5String seq;
    
    seqan::SequenceStream seqStream(argv[1]);
    readRecord(id, seq, seqStream);
    std::cout << id << '\t' << seq << '\n';
    
    return 0;
}