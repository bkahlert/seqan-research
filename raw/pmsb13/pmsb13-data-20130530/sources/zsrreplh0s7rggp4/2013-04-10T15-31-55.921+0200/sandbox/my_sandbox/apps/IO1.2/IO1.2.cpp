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

    seqan::CharString id1 = "seq1";
    seqan::Dna5String seq1 = "CGACGAACCGTT";

    if (writeRecord(seqStream, id1, seq1) != 0)
    {
        std::cerr << "ERROR: Could not write to file!\n";
        return 1;
    }

    seqan::CharString id2 = "seq2";
    seqan::Dna5String seq2 = "CGGTACCAGGTTCAAGTT";

    if (writeRecord(seqStream, id1, seq2) != 0)
    {
        std::cerr << "ERROR: Could not write to file!\n";
        return 1;
    } 
    
    
    return 0;
}