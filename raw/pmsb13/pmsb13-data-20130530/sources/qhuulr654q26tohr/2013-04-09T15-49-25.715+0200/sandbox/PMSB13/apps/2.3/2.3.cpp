#include <iostream>
#include <seqan/sequence.h>
#include <seqan/file.h>

using namespace seqan;

int main()
{
    StringSet<String<Dna5> > nucleotideSets;
    String<Dna5> seq1 ="ATATANGCGT";
    String<Dna5> seq2 ="AAGCATGANT";
    String<Dna5> seq3 ="TGAAANTGAC";
    appendValue(nucleotideSets,seq1);
    appendValue(nucleotideSets,seq2);
    appendValue(nucleotideSets,seq3);
    String<Dna5> refSeq="GATGCATGAT";
    StringSet<String<Dna5> > lesser;
    StringSet<String<Dna5> > greater;
    
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