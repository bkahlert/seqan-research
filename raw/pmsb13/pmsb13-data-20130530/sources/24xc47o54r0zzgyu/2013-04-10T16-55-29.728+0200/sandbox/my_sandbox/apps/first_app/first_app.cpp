#include <iostream>
#include <seqan/seq_io.h>

int main(int argc, char const ** argv)
{
    if (argc != 5)
    {
        std::cerr << "USAGE: INPUT.fa ID beginPos endPos";
        return 1;
    }

   	seqan::FaiIndex faiIndex;
	if (build(faiIndex, argv[1]) != 0)
	{
		std::cerr << "ERROR: Could not build the index!\n";
		return 1;
	}
	
	unsigned idx = 0;
    if (!getIdByName(faiIndex, argv[2], idx))
    {
        std::cerr << "ERROR: Unkown sequence-id: " << argv[2] << "\n";
        return 1;
    }

	unsigned beginPos;
	unsigned endPos;
	seqan::lexicalCast2(beginPos, argv[3]);
	seqan::lexicalCast2(endPos, argv[4]);
	if (beginPos>endPos)
	{
		std::cout << "Begin position must not be greater than end position" << std::endl;
	}

	seqan::Dna5String seq;
	if (readRegion(seq, faiIndex, idx, beginPos, endPos) != 0)
		std::cerr << "ERROR: Could not load" << argv[2]<<std::endl;
	std::cout << "Sequence ID: " << seq << std::endl;
	std::cout << "Position " << beginPos << " to " << endPos << std::endl;
	std::cout << seq << std::endl;

    return 0;
}