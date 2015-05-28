#include <iostream>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

int main(int argc, char const ** argv)
{
    if (argc >= 2)
    {
        seqan::CharString id;
        seqan::Dna5String seq;

        seqan::SequenceStream seqStream(argv[1]);
        readRecord(id, seq, seqStream);
        std::cout << id << '\t' << seq << '\n';

        return 0;
    }
    else
    {
        std::cout << "alles falsch";
        return 1;
    }
}
