#include <iostream>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

int main(int argc,const char * argv[])
{
    std::cout<<argv[0]<<std::endl;
    seqan::CharString id;
    seqan::Dna5String seq;

    seqan::SequenceStream seqStream(argv[1]);
    readRecord(id, seq, seqStream);
    std::cout << id << '\t' << seq << std::endl;

    return 0;
}
