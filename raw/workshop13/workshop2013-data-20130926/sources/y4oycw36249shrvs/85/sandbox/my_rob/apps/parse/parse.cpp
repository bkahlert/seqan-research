#include <iostream>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

int main(int argc, char const ** argv){
	
	if (argc < 2){
		std::cerr << "Error"<< std::endl;
		return 1;
	}
    seqan::CharString id;
    seqan::Dna5String seq;
    seqan::SequenceStream seqStream("exemple.fa");
    readRecord(id, seq, seqStream);
    std::cout << id << '\t' << seq << '\n';
    return 0;
}