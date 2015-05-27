#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <iostream>

using namespace seqan;


// Program entry point
int main()
	{
		// Each, a string of characters, of nucleotides and of amino acids.
		String<char> myText;
		String<Dna> myGenome;
		String<AminoAcid> myProtein;
		
		// More exotic: A string of character strings!
		String<String<char> > myStringList;
		
		String<char> test = "this is ";
		test += "a test.";
		std::cout << test << std::endl;
		std::cout << length(tes) << std::endl;
		
		return 0;
	}