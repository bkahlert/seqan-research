#include <iostream>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

int main(int argc, char**argv[])
{
	//if(argc<2){std::cerr << "nichts angegeben" << std::endl; return 1;}
	//seqan::String<char> path = "E:\example.fa";
    seqan::CharString id;
    seqan::Dna5String seq;
	//se
	seqan::SequenceStream seqStream("E:\example.fa");
	
	if (!isGood(seqStream))
    {
        std::cerr << "ERROR: Could not open the file.\n";
        return 1;
    }
	while(!atEnd(seqStream)){
		if (readRecord(id, seq, seqStream) != 0)
		{
			std::cerr << "ERROR: Could not read from example.fa!\n";
			return 1;
		}
		std::cout << id << '\t' << seq << '\n';
	}

	seqan::CharString id2;
	seqan::Dna5String seq2;
	seqan::CharString qly;
	seqan::SequenceStream seqStream2("E:\example.fq");
	if (!isGood(seqStream2))
    {
        std::cerr << "ERROR: Could not open the file.\n";
        return 1;
    }
	while(!atEnd(seqStream2)){
		if (readRecord(id2, seq2, qly, seqStream2) != 0)
		{
			std::cerr << "ERROR: Could not read from example.fq!\n";
			return 1;
		}
		std::cout << id2 << '\t' << seq2 << '\t' << qly << '\n';
	}
	// Write part
	seqan::SequenceStream seqStream3("E:\out.fq", seqan::SequenceStream::WRITE);
    if (!isGood(seqStream3))
    {
        std::cerr << "ERROR: Could not open the file.\n";
        return 1;
    }

    seqan::CharString id3 = "seq2";
    seqan::Dna5String seq3 = "CGAT";
	seqan::CharString qly3 = "IIII";

    if (writeRecord(seqStream3, id3, seq3, qly3) != 0)
    {
        std::cerr << "ERROR: Could not write to file!\n";
        return 1;
    }
    return 0;
}