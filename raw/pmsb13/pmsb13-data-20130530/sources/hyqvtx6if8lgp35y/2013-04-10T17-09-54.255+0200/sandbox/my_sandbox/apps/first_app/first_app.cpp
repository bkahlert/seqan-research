#include <iostream>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/stream.h>

int main(int argc, char const ** argv)
{
    if (argc != 5)
    {
        std::cerr << "Eingabefehler" << std::endl;
        return 1;
    }
    seqan::FaiIndex faiIndex;
    if (build(faiIndex, argv[1]) != 0)
    {
        std::cerr << "ERROR: Could not build FAI index for file." << std::endl;
        return 1;
    }
    seqan::CharString pathToFile = argv[1];
    append(pathToFile, ".fai");
    if (write(faiIndex, toCString(pathToFile)) != 0)
    {
        std::cerr << "ERROR: Could not write the index to file!" <<std::endl;
        return 1;
    }
    std::cout << "Datei " << pathToFile << " wurde erfolgreich erstellt" << std::endl;
	unsigned idx = 0;
    if (!getIdByName(faiIndex, argv[2], idx))
    {
        std::cerr << "ERROR: Index does not know about sequence " << argv[2] << std::endl;
        return 1;
    }
	
	unsigned begin;
	unsigned end;
	if (!seqan::lexicalCast2(begin, argv[3]))
    {
        std::cerr << "ERROR " << argv[3] << std::endl;
        return 1;
    }
    if (!seqan::lexicalCast2(end, argv[4]))
    {
        std::cerr << "ERROR" << argv[4] << std::endl;
        return 1;
    }	
	if(begin > sequenceLength(faiIndex, argv[2]))	
		begin = sequenceLength(faiIndex, argv[2]);	
	if(begin > end)
		end = begin;
	if(end > sequenceLength(faiIndex, argv[2]))
		end = sequenceLength(faiIndex, argv[2]);
    if (readRegion(infix, faiIndex, idx, begin, end) != 0)
    {
        std::cerr << "ERROR" << std::endl;
        return 1;
    }
    std::cout << infix << std::endl;
    return 0;
}
