#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/align.h>
#include <iostream>
using namespace seqan;


int main(int, char **) {

	typedef String<char> TSequence;  // sequence type
    typedef Align<TSequence,ArrayGaps>  TAlign;     // align type
	typedef Row<TAlign>::Type TRow;  

	String<char> str1 = "blablablu";
	String<char> str2 = "abab";

	TAlign align;
    resize(rows(align), 2);
    assignSource(row(align,0),str1);
    assignSource(row(align,1),str2);

	int score = globalAlignment(align,Score<int>(EditDistance),TBottom);
	
	::std::cout << align << std::endl;
	::std::cout << "Score: " << score << std::endl;

    return 1;
}
