#include <iostream>
#include <seqan/sequence.h>
#include <seqan/file.h>

using namespace seqan;

void first_scan()
{
    String<Dna5> nucleotides = "AGTCGTGNNANCT";
	String<Dna5> g="G";
	String<Dna5> lesser;
	String<Dna5> greater;
    
    for (unsigned i = 0; i < length(nucleotides); ++i){
		if(nucleotides[i] < g)
			appendValue(lesser, nucleotides[i]);
		else
			if(nucleotides[i] != g)
				appendValue(greater, nucleotides[i]);
    }
    std::cout << "All nucleotides: " << nucleotides << "\n" << "Unselected nucleotides: " << g << "\n" << "Lesser nucleotides: "<< lesser << "\n" << "Greater nucleotides: " << greater << std::endl;
}

int main()
{

	String<Dna5String> nucleotides;
	append(nucleotides,"ATATANGCGT");
	append(nucleotides,"AAGCATGANT");
	append(nucleotides,"TGAAANTGAC");

	String<Dna5String> g="GATGCATGAT";
	String<Dna5String> lesser;
	String<Dna5String> greater;
    
    for (unsigned i = 0; i < length(nucleotides); ++i){
		if(nucleotides[i] < g)
			appendValue(lesser, nucleotides[i]);
		else
			if(nucleotides[i] != g)
				appendValue(greater, nucleotides[i]);
    }

	 std::cout << "All nucleotides: " << nucleotides << "\n" << "Unselected nucleotides: " << g << "\n" << "Lesser nucleotides: "<< lesser << "\n" << "Greater nucleotides: " << greater << std::endl;



	return 0;
}