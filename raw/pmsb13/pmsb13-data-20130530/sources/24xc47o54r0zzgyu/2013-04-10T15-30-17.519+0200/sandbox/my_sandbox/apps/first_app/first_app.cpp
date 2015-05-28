#include <iostream>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

int main(int argc, char const ** argv)
{
    if (argc < 2)
    {
        std::cerr << "USAGE: basic_seq_io_example FILENAME\n";
        return 1;
    }

    seqan::SequenceStream seqStream(argv[1], seqan::SequenceStream::WRITE);
    if (!isGood(seqStream))
    {
        std::cerr << "ERROR: Could not open the file.\n";
        return 1;
    }

  	seqan::StringSet<seqan::Dna5String> seqs;
	seqan::StringSet<seqan::CharString> ids;
	seqan::appendValue(ids, "seq1");
	seqan::appendValue(ids, "seq2");
	seqan::appendValue(seqs, "seq2");
	seqan::appendValue(seqs, "GATACA");

	for(unsigned i=0;i<2;++i)
	{
		if (writeRecord(seqStream, ids[i], seqs[i]) != 0)
		{
			std::cerr << "ERROR: Could not write to file!\n";
			return 1;
		}
	}
    return 0;
}