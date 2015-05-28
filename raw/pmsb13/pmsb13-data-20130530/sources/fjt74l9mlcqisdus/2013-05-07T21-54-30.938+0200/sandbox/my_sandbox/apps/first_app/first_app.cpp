// BLASTX IMPLEMENTIERUNG VON ANNKATRIN BRESSIN UND MARJAN FAIZI
// SOFTWAREPROJEKT VOM 2.4. - 29.5.2012
// VERSION VOM 

#include <iostream>
#include <seqan/arg_parse.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/find.h>
#include <seqan/index.h>
#include <seqan/basic.h>
#include <cstring>
#include <seqan/index_extras.h>
#include <seqan/score.h>
#include <algorithm>
#include <typeinfo>
#include <seqan/align.h>

using namespace std;
using namespace seqan;

int main(){

typedef String<char> TSequence;                 // sequence type
typedef Align<TSequence, ArrayGaps> TAlign;     // align type

TSequence seq1 = "GARFIELDTHECAT";
TSequence seq2 = "GARFIELD";

TAlign align;
resize(rows(align),2);
assignSource(row(align,0),seq1);
assignSource(row(align,1),seq2);
int score = globalAlignment(align,Score<int,Simple>(1,-1,-1), AlignConfig<true, true, true, true>(), NeedlemanWunsch());

cout << "Score: " << score <<endl;
cout << align <<endl;
cout<< row(align, 0)[0] <<endl;



    return 0;
}
