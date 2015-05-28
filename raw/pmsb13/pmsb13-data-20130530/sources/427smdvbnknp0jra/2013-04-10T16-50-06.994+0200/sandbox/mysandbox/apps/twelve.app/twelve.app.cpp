#include <iostream>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

int main(int argc,char argv)
{
	/*if(argc<2)
	{	
		std::cout << "-ERROR-" << std::endl;
		return 1;
	}*/

    seqan::CharString id;
    seqan::Dna5String seq;

    seqan::SequenceStream seqStream("D:\Informatik\Development\example.fa");
	if(!isGood(seqStream))
	{
		std::cout << "-ERROR- Could not open" << std::endl;
		return 1;
	}
	readRecord(id, seq, seqStream);
	if(readRecord(id, seq, seqStream)!=0)
	{
		std::cout << "-ERROR- Could not read" << std::endl;
		return 1;
	}
	
    std::cout << id << '\t' << seq << '\n';

    return 0;
}
