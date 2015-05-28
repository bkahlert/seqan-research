#include <iostream>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

int main(int argc, char * argv[])
{
   String <char> file;
   if(argc>0){
   file = argv[1];
   }
    std::cout<<argv[0]<<std::endl;
    seqan::CharString id;
    seqan::Dna5String seq;

    seqan::SequenceStream seqStream(file);
    readRecord(id, seq, seqStream);
    std::cout << id << '\t' << seq << std::endl;

    return 0;
}
