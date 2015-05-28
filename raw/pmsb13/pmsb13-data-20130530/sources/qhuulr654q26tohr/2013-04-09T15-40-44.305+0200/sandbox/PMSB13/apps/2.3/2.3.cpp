#include <iostream>
#include <seqan/sequence.h>
#include <seqan/file.h>

using namespace seqan;

int main()
{
    StringSet<Dna5> nucleotideSets;
    appendValue(nucleotideSets,new Dna5("ATATANGCGT"));
    appendValue(nucleotideSets,"AAGCATGANT");
    appendValue(nucleotideSets,"TGAAANTGAC");
    String<Dna5> refSeq="GATGCATGAT";
    StringSet<Dna5> lesser;
    StringSet<Dna5> greater;
    
    for (unsigned i = 0; i < length(nucleotideSets); ++i){
	Lexical<> comp(nucleotideSets[i],refSeq);
	if (isLess(comp))
	    appendValue(lesser, nucleotideSets[i]);
	if (isGreater(comp))
	    appendValue(greater, nucleotideSets[i]);
    }
    //std::cout << "lesser: " << lesser << std::endl;
    //std::cout << "greater: " << greater << std::endl;
    return 0;
}