#include <iostream>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

template <typename Tstring>
void build_fai(TString const & pathToFile)
{
	seqan::FaiIndex faiIndex;
	int res = build(faiIndex, pathToFile);
	if (res != 0)
    	std::cerr << "ERROR: Could not build the index!\n";
}

int main(int argc, char const ** argv)
{
	if(argc > 2)
	{
		std::cerr << "Eingabefehler: Zu viele Argumente" << std::endl;		
		return 1;
	}	
	build_fai(arg[1]);	
	int res = write(faiIndex, argv[1]);
	if (res != 0)
    	std::cerr << "ERROR: Could not write the index to file!\n";

	std::cout << "Datei wurde erfolgreich erstellt" << std::endl;
    return 0;
}
