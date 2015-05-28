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

	if (!getIdByName(faiIndex, argv[2], argv[3]))
    std::cerr << "ERROR: FAI index has no entry for chr1.\n";
	
    return 0;
}