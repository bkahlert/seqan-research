// Program part for parsing the input
#include <iostream>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

using namespace seqan;

int main(int argc, char const ** argv)
{
	if (argc < 2)
	{
		std::cerr << "Please input the path to the input file as first argument.\n";
		return 1;
	}
	
	StringSet<seqan::CharString> ids;
	StringSet<seqan::Dna5String> seqs;
	StringSet<seqan::CharString> quals;
	int res = 0;
	SequenceStream seqStream(argv[1]);
	unsigned reads = 100;									//number of records to be read
	resize(ids, reads);
	resize(seqs, reads);
	resize(quals, reads);

	if (!isGood(seqStream))
	{
		std::cerr << "Error while opening the input file.\n";
		return 1;
	}

	res = readBatch(ids, seqs, quals, seqStream, reads);
	if(res != 0)		
	{
		std::cerr << "Error while reading the sequences.\n";
		return 1;
	}

	for (unsigned i = 0; i < reads; ++i)
		std::cout << ids[i] << '\t' << seqs[i] << '\t' << quals[i] << '\n'; //Testausgabe

	return 0;
}