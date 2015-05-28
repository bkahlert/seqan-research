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

	seqan::StringSet<seqan::CharString> ids;
	seqan::StringSet<seqan::Dna5String> seqs;
	appendValue(ids, "seq1");
	appendValue(ids, "seq2");
	appendValue(seqs, "CGAT");
	appendValue(seqs, "CATAG");

	if (writeAll(seqStream, ids, seqs) != 0)
    {
        std::cerr << "ERROR: Could not write to file!\n";
        return 1;
    }

    return 0;
}
