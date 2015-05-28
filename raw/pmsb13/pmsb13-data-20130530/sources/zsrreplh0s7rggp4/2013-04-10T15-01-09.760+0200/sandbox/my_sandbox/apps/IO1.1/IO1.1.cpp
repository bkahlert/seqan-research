#include <iostream>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

int main(int argc,char const ** argv)
{
   
    if (argc < 2)
    {
        std::cerr << "Error\n";
        return 1;
    }
    
    seqan::CharString id;
    seqan::Dna5String seq;
    for (int i=1;i<length(seqStream);i++){
    seqan::SequenceStream seqStream(argv[i]);
    while(!atEnd(SequenceStream)){readRecord(id, seq, seqStream);
    std::cout << id << '\t' << seq << '\n';}

    if (!isGood(seqStream))
    {
        std::cerr << "ERROR: Could not open the file.\n";
        return 1;
    }
    if (readRecord(id, seq, seqStream) != 0)
    {
        std::cerr << "ERROR: Could not read from example.fa!\n";
        return 1;
    }
    
    return 0;
}