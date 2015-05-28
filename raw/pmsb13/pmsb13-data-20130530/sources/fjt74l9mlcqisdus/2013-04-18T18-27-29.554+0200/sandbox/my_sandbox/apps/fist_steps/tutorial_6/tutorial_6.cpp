#include <iostream>
#include <seqan/sequence.h>
#include <seqan/file.h>

using namespace seqan;

int main()
{
    String<String<Dna5>> listOfSeq;
	append(listOfSeq,"ATATANGCGT");
	append(listOfSeq,"AAGCATGANT");
	append(listOfSeq,"TGAAANTGAC");
    String<Dna5> referenz="GATGCATGAT";
    
	for (unsigned i = 0; i < length(nucleotides); ++i){
        Lexical<> comp(listOfSeq[i], referenz);
		appendValue(selected, nucleotides[i]);
    }
    std::cout << "Selected nucleotides: " << selected << std::endl;
    system("Pause");
	return 0;
}
