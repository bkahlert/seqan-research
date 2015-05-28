#include <iostream>
#include <seqan/sequence.h>
#include <seqan/file.h>

using namespace seqan;

int main()
{
    String<Dna5String> dnas = (Dna5String) "AGTCGTGNNANCT";
	dnas += "AAGCATGANT";
	dnas += "TGAAANTGAC";
	Dna5String reference = "GATGCATGAT";

	std::cout << dnas; 

    String<Dna5String> lesser;
	String<Dna5String> greater;

	for (unsigned i=0;0<length(dnas)-1;++i)
	{
	Lexical<> comp (dnas[i], reference);
	if (isLess(comp))
		lesser += dnas[i];
	else if (isGreater(comp))
		greater += dnas[i];
	}

    std::cout << "Lesser Sequences: " << lesser << std::endl;
	std::cout << "Greater Sequences: " << greater << std:: endl;
    return 0;
}