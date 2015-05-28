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
    CharString qua;
    seqan::SequenceStream seqStream(argv[1], seqan::SequenceStream::WRITE);

    if (!isGood(seqStream))
    {
        std::cerr << "ERROR: Could not open the file.\n";
        return 1;
    }

//    while (!atEnd(seqStream))
//    {
//        if (readRecord(id, seq, qua, seqStream) != 0)
//        {
//            std::cerr << "ERROR: Could not read from example.fa!\n";
//            return 1;
//        }
//        std::cout << id << '\t' << seq << '\t'<< qua <<'\n';
//    }
    seqan::CharString id = "seq1";
    seqan::Dna5String seq = "CGAT";
    if (readRecord(id, seq, qua, seqStream) != 0)
            {
                std::cerr << "ERROR: Could not write to example.fa!\n";
                return 1;
            }
    seqan::CharString id = "seq2";
    seqan::Dna5String seq = "AAAAAA";
    if (readRecord(id, seq, qua, seqStream) != 0)
            {
                std::cerr << "ERROR: Could not write to example.fa!\n";
                return 1;
            }

    return 0;
}

