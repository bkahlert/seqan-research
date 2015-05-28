#include <iostream>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

int main()
{
    seqan::CharString id;
    seqan::Dna5String seq;

    seqan::SequenceStream seqStream("/home/stefan/Dokumente/subversion/seqan-trunk/sandbox/my_sandbox/apps/bsio_a1/example.fa");
    readRecord(id, seq, seqStream);
    std::cout << id << '\t' << seq << '\n';

    return 0;
}

