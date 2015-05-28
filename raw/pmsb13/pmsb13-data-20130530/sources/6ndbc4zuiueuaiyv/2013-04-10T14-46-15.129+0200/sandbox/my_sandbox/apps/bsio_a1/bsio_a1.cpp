#include <iostream>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <string>
using namespace std;

int main()
{
    if (argc < 2) {
        cout << "Usage: " << argv[0] << " <fasta-file>" << endl;
        return 1
    }
    seqan::CharString id;
    seqan::Dna5String seq;

    seqan::SequenceStream seqStream(argv[1]);
    readRecord(id, seq, seqStream);
    std::cout << id << '\t' << seq << '\n';

    return 0;
}

