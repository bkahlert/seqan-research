#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <iostream>
#include <seqan/align.h>

using namespace seqan;

//template <typename TAlphabet>
typedef String<char> TSequence;                 // sequence type
template <typename TSequence>

void align_sequences(TSequence & seqa, TSequence & seqb)
{
//	TSequence seqa = "aa";
//	TSequence seqb = "ab";
	typedef Align<TSequence,ArrayGaps> TAlign;      // align type
    //typedef Row<TAlign>::Type TRow;     
	
	TAlign align;
    //resize(rows(align), 2);
    //assignSource(row(align,0),seqa);
    //assignSource(row(align,1),seqb);
	//::std::cout << align;	
}

int main()
{

align_sequences("aaa", "baa")
}