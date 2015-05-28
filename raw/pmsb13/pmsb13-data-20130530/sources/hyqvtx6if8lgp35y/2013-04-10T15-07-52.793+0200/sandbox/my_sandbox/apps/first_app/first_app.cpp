#include <iostream>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

int main(int argc, char const ** argv)
{
	if(argc > 2)
	{
		std::cerr << "Eingabefehler: Zu viele Argumente" << std::endl;
		return 1;
	}    
	seqan::CharString id;
    seqan::Dna5String seq;
	seqan::CharString qual;
    seqan::SequenceStream seqStream(argv[1]);
	if(!isGood(seqStream))
	{
		std::cerr << "Error:Could not open file." << std::endl;
		return 1;
	}    
	while(!atEnd(seqStream))
	{	
		if(readRecord(id, seq, qual; seqStream) != 0)
		{
			std::cerr << "Error:Could not open file." << std::endl;
			return 1;
		}
    
		std::cout << id << '\t' << seq << '\t' << qual << '\n';
	}    

	return 0;
}
