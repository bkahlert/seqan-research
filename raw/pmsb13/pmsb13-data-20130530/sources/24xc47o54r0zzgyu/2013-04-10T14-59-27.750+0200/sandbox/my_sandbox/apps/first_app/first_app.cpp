#include <iostream>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

int main(int argc, char* argv[])
{	if (argc < 2)
	{
		std::cout << "Kein Parameter uebergeben" <<std::endl;
		return 1;
	}
	else
	{
		seqan::SequenceStream seqStream(argv[1]);
		if (!isGood(seqStream))
		{
			std::cerr << "ERROR: Could not open the file.\n";
			return 1;
		}
		
		for(;!atEnd(seqStream);)
		{
			seqan::CharString id;
			seqan::Dna5String seq;
			if (readRecord(id, seq, seqStream) != 0)
			{
			std::cerr << "ERROR: Could not read from example.fa!\n";
			return 1;
			}
			readRecord(id, seq, seqStream);
			std::cout << id << '\t' << seq << '\n';
		}

	}
    return 0;
}