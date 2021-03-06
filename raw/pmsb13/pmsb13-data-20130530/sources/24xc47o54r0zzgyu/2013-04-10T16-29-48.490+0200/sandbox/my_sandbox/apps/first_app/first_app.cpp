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
	
	unsigned idx = (unsigned int) argv[3];
	if (!getIdByName(faiIndex, argv[2], idx))
		std::cerr << "ERROR: FAI index has no entry matching " << argv[2] <<std::endl;
	
	seqan::Dna5String seq;
	unsigned idxend = (unsigned int) argv[3];
	if (readRegion(seq, faiIndex, idx, 0, idxend-idx) != 0)
		std::cerr << "ERROR: Could not load chr1.\n";
	std::cout << seq << std::endl;

    return 0;
}