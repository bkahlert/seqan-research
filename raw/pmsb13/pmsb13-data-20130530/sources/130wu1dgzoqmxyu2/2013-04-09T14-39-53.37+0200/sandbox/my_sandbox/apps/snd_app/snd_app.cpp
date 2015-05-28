#include <iostream>
#include <seqan/sequence.h>
#include <seqan/file.h>

using namespace seqan;

void biglessThan(String<Dna5> bla, String<Dna5> pattern){
	String<Dna5> greaterThan;
	String<Dna5> lesserThan;
	for (unsigned i = 0; i < length(pattern); ++i){
		if (pattern[i] > bla[i]){
			appendValue(lesserThan, bla[i]);
		}
		else {
			appendValue(greaterThan, bla[i]);
		}
	}
	std::cout<<"greater: "<<greaterThan<<std::endl;
	std::cout<<"lesser: "<<lesserThan<<std::endl;
}

int main()
{
    String<Dna5> nucleotides = "AGTCGTGNNANCT";
    String<Dna5> lesser;
    String<Dna5> greater;
    String<Dna5> gama = "ATATANGCGT";
    String<Dna5> beta = "AAGCATGANT";
    String<Dna5> alpha = "TGAAANTGAC";
    String<Dna5> pattern = "GATGCATGAT";

    biglessThan(alpha, pattern);
    biglessThan(beta, pattern);
    biglessThan(gama, pattern);

}
