#include <iostream>
#include <fstream>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>  // CharString, ...
#include <seqan/file.h> // to stream a CharString into cout

using namespace seqan;

bool readFASTA(char const * path, CharString  &id, Dna5String &seq){

	std::fstream in(path, std::ios::binary | std::ios::in);
	RecordReader<std::fstream, SinglePass<> > reader(in);

	if (readRecord(id, seq, reader, seqan::Fasta()) == 0){
		return true;
	}
	else{
		return false;
	}
}

int main(int argc, char const ** argv){


	CharString idOne;
	Dna5String seqOne;
	CharString idTwo;
	Dna5String seqTwo;

	if (!readFASTA(argv[1], idOne, seqOne)){
		return 1;
	}

	if (!readFASTA(argv[2], idTwo, seqTwo)){
		return 1;
	}

	std::cout << idOne << "\n" << seqOne;
	return 0;
}
