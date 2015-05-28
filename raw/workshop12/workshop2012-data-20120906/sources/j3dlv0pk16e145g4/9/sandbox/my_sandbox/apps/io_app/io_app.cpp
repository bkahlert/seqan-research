#include <iostream>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

int main(int argc, char const ** argv)
{
    seqan::CharString id;
    seqan::Dna5String seq;
    if (argc < 2)
    {
        std::cerr << "Please give an input file." << std::endl;
        return 1;
    }
    seqan::SequenceStream seqStream(argv[1]);
    readRecord(id, seq, seqStream);
    std::cout << id << '\t' << seq << '\n';
    return 0;
}