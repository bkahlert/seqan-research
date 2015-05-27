#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <seqan/find_motif.h>
#include <seqan/file.h>
#include <iostream>
#include <seqan/align.h>

using namespace seqan;

int main (){

	typedef String<Dna> TDnaSeq;
	typedef Align<TDnaSeq, ArrayGaps> TAlign;


	TDnaSeq seq1 = "acgtacgtact";
	TDnaSeq seq2 = "actactacgt";
	TAlign align;
	resize(rows(align),2);
	assignSource(row(align,0),seq1);
	assignSource(row(align,1),seq2);

	::std::cout<<align;

	int score = globalAlignment(align,Score<int>(1,-1,-1,-1), Hirschberg());

	::std::cout<<align;



	return 0;
}