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
	resize(nucleotides,4);
	nucleotides[0]="ATATANGCGT";
	nucleotides[1]="AAGCATGANT";
	nucleotides[2]="TGAAANTGAC";

	String<Dna5String> g="GATGCATGAT";
	String<int> lesser;
	String<int> greater;
    
    for (unsigned i = 0; i < length(nucleotides); ++i){
		if(nucleotides[i] < g)
			appendValue(lesser, i);
		else
			if(nucleotides[i] != g)
				appendValue(greater, i);
    }

	std::cout << "All nucleotides: " ;

	for (unsigned i = 0; i < length(nucleotides); ++i){ 
		std::cout << nucleotides[i] << " ";
	}
	std::cout << "\nUnselected nucleotides: " << g << "\n" << "Lesser nucleotides: ";
	for (unsigned i = 0; i < length(lesser); ++i){
		std::cout << nucleotides[lesser[i]] << " " ;
	}
	std::cout << "\nGreater nucleotides: " ;
	for (unsigned i = 0; i < length(greater); ++i){
		std::cout << nucleotides[greater[i]] << " " ;
	}
	



	return 0;
}