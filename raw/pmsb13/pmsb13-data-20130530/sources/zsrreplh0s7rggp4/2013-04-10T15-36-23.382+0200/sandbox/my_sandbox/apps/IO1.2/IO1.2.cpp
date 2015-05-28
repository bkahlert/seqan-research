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

    seqan::CharString id2 = "seq2";
    seqan::Dna5String seq2 = "CGGTACCAGGTTCAAGTT";
    
    seqan::StringSet<seqan::CharString> ids;
    seqan::StringSet<seqan::Dna5String> seqs;
    seqan::StringSet<seqan::CharString> quals;
    
    seqs[0]=seq1;
    seqs[1]=seq2;
    ids[0]=id1;
    ids[1]=id2;
    
    if (writeAll(seqStream, ids, seqs) != 0)
    {
        std::cerr << "ERROR: Could not write to file!\n";
        return 1;
    }

  
    
    return 0;
}