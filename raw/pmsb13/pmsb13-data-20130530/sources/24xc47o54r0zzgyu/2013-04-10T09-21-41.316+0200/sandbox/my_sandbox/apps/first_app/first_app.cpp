#include <iostream>
#include <seqan/align.h>

using namespace seqan;

typedef String<Dna> TSequence;
typedef Align<TSequence,ArrayGaps> TAlign;
typedef Row<TAlign>::Type TRow;

TSequence seq1 = "ACGTCACCTC";
TSequence seq2 = "ACGGGCCTATC";

TAlign align;
resize(rows(align), 2);
assignSource(row(align,0),seq1);
assignSource(row(align,1),seq2);


int main()
{
   
    return 0;
}