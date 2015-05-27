#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <iostream>
#include <seqan/align.h>

using namespace seqan;

//template <typename TAlphabet>
//template <typename TSequence>
void align_sequences()
//void align_sequences(TSequence & seqa, TSequence & seqb)
{
	typedef String<char> TSequence;                 // sequence type

	typedef Align<TSequence, ArrayGaps> TAlign;      // align type
    typedef Row<TAlign>::Type TRow;     
	
	TAlign align;
    
	//	TSequence seqa = "aa";
	//	TSequence seqb = "ab";

	//resize(rows(align), 2);
    //assignSource(row(align,0),seqa);
    //assignSource(row(align,1),seqb);
	//::std::cout << align;	
}

int main()
{
	align_sequences();
//	align_sequences("aaa", "baa");
}