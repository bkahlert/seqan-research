#include <iostream>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

int main(int argc, char const ** argv)
{
    if (argc > 2)
    {
        std::cerr << "Eingabefehler: Zu viele Argumente" << std::endl;
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
    if (write(faiIndex, toCString(argv[1]), pathToFile) != 0)
    {
        std::cerr << "ERROR: Could not write the index to file!" <<std::endl;
        return 1;
    }
    std::cout << "Datei " << pathToFile << " wurde erfolgreich erstellt" << std::endl;
    return 0;
}
