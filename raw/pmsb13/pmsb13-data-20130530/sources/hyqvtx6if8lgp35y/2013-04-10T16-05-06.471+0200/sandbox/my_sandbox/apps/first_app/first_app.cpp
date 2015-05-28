#include <iostream>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

template <typename TString>
seqan::FaiIndex build_fai(TString const & pathToFile)
{
	seqan::FaiIndex faiIndex;
	int res = build(faiIndex, pathToFile);
	if (res != 0)
    	std::cerr << "ERROR: Could not build the index!\n";
	return faiIndex;
}

int main(int argc, char const ** argv)
{
	if(argc > 2)
	{
		std::cerr << "Eingabefehler: Zu viele Argumente" << std::endl;		
		return 1;
	}	
	seqann::FaiIndex index = build_fai(argv[1]);	
	int res = write(index, argv[1]);
	if (res != 0)
    	std::cerr << "ERROR: Could not write the index to file!\n";

	std::cout << "Datei wurde erfolgreich erstellt" << std::endl;
    return 0;
}
