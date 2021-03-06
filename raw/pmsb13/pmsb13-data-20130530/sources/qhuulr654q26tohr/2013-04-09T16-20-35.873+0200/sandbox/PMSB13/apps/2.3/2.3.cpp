#include <iostream>
#include <seqan/sequence.h>
#include <seqan/file.h>

using namespace seqan;
using namespace std;


void print(StringSet<String<Dna5> > a){
    for (int i=0;i<length(a);++i){
	cout << a[i];
	cout << ", ";}
	cout << endl;
}


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
    //std::cout << "lesser: " << std::endl;
    //print(lesser);
    //std::cout << "greater: " << std::endl;
    cout << length(greater) << endl;
    //print(greater);
    return 0;
}