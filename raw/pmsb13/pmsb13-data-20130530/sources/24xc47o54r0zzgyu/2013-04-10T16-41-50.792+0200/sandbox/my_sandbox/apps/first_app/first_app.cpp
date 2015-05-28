#include <iostream>
#include <seqan/seq_io.h>

int main(int argc, char const ** argv)
{
    if (argc != 5)
    {
        std::cerr << "USAGE: INPUT.fa ID Begin End";
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

	unsigned beginpos;
	unsigned endpos;
	seqan::lexicalCast2(beginpos, argv[3]);
	seqan::lexicalCast2(endpos, argv[4]);

	seqan::Dna5String seq;
	if (readRegion(seq, faiIndex, idx, 0, endpos-beginpos) != 0)
		std::cerr << "ERROR: Could not load" << argv[2]<<std::endl;
	std::cout << seq << std::endl;

    return 0;
}