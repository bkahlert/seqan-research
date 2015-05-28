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

    seqan::CharString id;
    seqan::Dna5String seq;

    seqan::SequenceStream seqStream(argv[1]);

    if (!isGood(seqStream))
    {
        std::cerr << "ERROR: Could not open the file.\n";
        return 1;
    }
    if (readRecord(id, seq, seqStream) != 0)
    {
        std::cerr << "ERROR: Could not read from example.fa!\n";
        return 1;
    }
    Iterator<CharString, Rooted>::Type itId = begin(id);
    Iterator<Dna5String, Rooted>::Type itSeq = begin(seq);
    for (; !atEnd(itId); itId++, itSeq++)
    {
        std::cout << getValue(itId) << '\t' << getValue(itSeq) << '\n';
    }

    return 0;
}

