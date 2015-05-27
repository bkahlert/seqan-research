#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/align.h>
#include <iostream>
using namespace seqan;


int main(int, char **) {

	typedef String<Dna> TSequence;  // sequence type
    typedef Align<TSequence,ArrayGaps>  TAlign;     // align type

	String<Dna> str1 = "acgtacgtact";
	String<Dna> str2 = "actactacgt";

	TAlign align;
    resize(rows(align), 2);
    assignSource(row(align,0),str1);
    assignSource(row(align,1),str2);

	int score = globalAlignment(align,Score<int>(1,-1,-1,-1));

	::std::cout << align;

	::std::cout << ::std::endl << "ViewPos1: ";
    for(unsigned i = 0; i < length(source(row1)); ++i){
		if(isGap(source(row1),i)){
			::std::cout << toViewPosition(row1, i) << ",";
		}
	}
    ::std::cout << ::std::endl << "ViewPos2: ";
    for(unsigned i = 0; i < length(source(row2)); ++i){
		if(isGap(source(row2),i)){
			::std::cout << toViewPosition(row2, i) << ",";
		}
	}
    ::std::cout << ::std::endl;

    return 1;
}
