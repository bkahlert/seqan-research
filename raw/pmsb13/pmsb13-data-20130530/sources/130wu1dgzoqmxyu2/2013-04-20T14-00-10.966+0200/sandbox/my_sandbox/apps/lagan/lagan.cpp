#include <iostream>
#include <fstream>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>  // CharString, ...
#include <seqan/file.h> // to stream a CharString into cout

using namespace seqan;

void readFASTA(char const * path, CharString * id, Dna5String * seq){

	std::fstream in(path, std::ios::binary | std::ios::in);
	RecordReader<std::fstream, SinglePass<> > reader(in);

	readRecord(id, seq, reader, seqan::Fasta());
	std::cout << id;
}

int main(int argc, char const ** argv){

	CharString idOne;
	Dna5String seqOne;
	std::cout<<argv[1];
	readFASTA(argv[1], idOne, seqOne);

	std::cout << idOne << "\n" << seqOne;
	return 0;
}
