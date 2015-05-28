#include <iostream>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

int main()
{
    if (argc < 2) {
        cout << "Usage: " << argv[0] << " <fasta-file>" << endl;
        return 1
    }
    seqan::CharString id;
    seqan::Dna5String seq;

    seqan::SequenceStream seqStream("/home/stefan/Dokumente/subversion/seqan-trunk/sandbox/my_sandbox/apps/bsio_a1/example.fa");
    readRecord(id, seq, seqStream);
    std::cout << id << '\t' << seq << '\n';

    return 0;
}

