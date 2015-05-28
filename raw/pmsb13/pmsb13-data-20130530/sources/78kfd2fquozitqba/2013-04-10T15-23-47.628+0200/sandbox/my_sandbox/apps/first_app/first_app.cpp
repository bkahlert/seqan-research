#include <iostream>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/basic.h>
#include <seqan/file.h>
using namespace seqan;
int main(int argc, char const ** argv)
{
    if (argc < 2)
    {
        std::cerr << "USAGE: basic_seq_io_example FILENAME\n";
        return 1;
    }

    seqan::StringSet<seqan::CharString> ids;
    seqan::StringSet<seqan::Dna5String> seqs;
    CharString qua;
    seqan::SequenceStream seqStream(argv[1], seqan::SequenceStream::WRITE);

    if (!isGood(seqStream))
    {
        std::cerr << "ERROR: Could not open the file.\n";
        return 1;
    }

    appendValue(ids, "seq 1");
    appendValue(ids, "seq 2");
    appendValue(seqs, "CGAT");
    appendValue(seqs, "AAAAAA");

    if (writeAll(seqStream, ids, seqs) != 0)
    {
        std::cerr << "ERROR: Could not write to example.fa!\n";
        return 1;
    }
    return 0;
}

