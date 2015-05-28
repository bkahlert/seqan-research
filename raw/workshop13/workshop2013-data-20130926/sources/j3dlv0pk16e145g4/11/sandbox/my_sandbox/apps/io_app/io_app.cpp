#include <iostream>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

int main(int argc, char const ** argv)
{
    seqan::CharString id;
    seqan::Dna5String seq;
    seqan::CharSting qual;
    if (argc < 2)
    {
        std::cerr << "Please give an input file." << std::endl;
        return 1;
    }
    seqan::SequenceStream seqStream(argv[1]);
    if (!isGood(seqStream))
    {
        std::cerr << "Could not open file " << argv[1] << std::endl;
        return 1;
    }
    while (!atEnd(seqStream))
    {
        if (readRecord(id, seq, qual, seqStream) != 0)
        {
            std::cerr << "Could not read record" << std::endl;
            return 1;
        }
        std::cout << id << '\t' << seq << '\n';
    }
    return 0;
}