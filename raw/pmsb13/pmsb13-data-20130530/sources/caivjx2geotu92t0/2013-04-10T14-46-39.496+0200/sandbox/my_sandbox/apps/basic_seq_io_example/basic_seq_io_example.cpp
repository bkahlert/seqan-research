#include <iostream>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

int main(int argc, char**argv[])
{
	//if(argc )
    seqan::CharString id;
    seqan::Dna5String seq;
	//se
	seqan::SequenceStream seqStream("C:\Daten\example.fa");
	readRecord(id, seq, seqStream);
	/*if (!isGood(seqStream))
    {
        std::cerr << "ERROR: Could not open the file.\n";
        return 1;
    }
    if (readRecord(id, seq, seqStream) != 0)
    {
        std::cerr << "ERROR: Could not read from example.fa!\n";
        return 1;
    }*/

    
    
    std::cout << id << '\t' << seq << '\n';

    return 0;
}