#include <iostream>
#include <seqan/align.h>

using namespace seqan;

void assignment1(){
    typedef String<char> TSequence;                 // sequence type
    typedef Align<TSequence,ArrayGaps> TAlign;      // align type
    typedef Row<TAlign>::Type TRow;                 // gapped sequence type

    TSequence seq1 = "ACGTCACCTC";
    TSequence seq2 = "ACGGGCCTATC";

    TAlign align;
    resize(rows(align), 2);
    assignSource(row(align,0),seq1);
    assignSource(row(align,1),seq2);

    std::cout << align;
    TRow &row1 = row(align,0);
    TRow &row2 = row(align,1);
    insertGap(row1,5);
    insertGaps(row1,2,2);
    insertGaps(row2,9,2);
    std::cout << align;

    typedef Iterator<TRow>::Type TRowIterator;
    TRowIterator it = begin(row1);
    TRowIterator itEnd = end(row1);

    unsigned gaps = 0;
    unsigned count;
    for (; it != itEnd;){
    	count = countGaps(it);
    	gaps += count;
    	it += count > 0 ? count : 1;
    }

    it = begin(row2);
    itEnd = end(row2);

    for (; it != itEnd;){
    	count = countGaps(it);
    	gaps += count;
    	it += count > 0 ? count : 1;
    }

    std::cout << gaps << std::endl;
}

void assignment2(){
}

// approximate pattern matching
void assignment5(){
    typedef String<char> TSequence;                // sequence type
    typedef Align<TSequence,ArrayGaps> TAlign;      // align type

    TSequence text = "MISSISSIPPIANDMISSOURI";
    TSequence pattern = "SISSI";

    // Iterieren über String und alignieren alle Positionen nach Edit-Distanz mit Muster
    for (unsigned i=0; i < length(text); ++i){
        TAlign align;
        resize(rows(align), 2);
        assignSource(row(align,0),infix(text, i, i+ length(pattern)));
        assignSource(row(align,1),pattern);

        int score = globalAlignment(align, Score<int,Simple>(0,-1,-1));

        ::std::cout << "Score: " << score << ::std::endl;
        ::std::cout << align << ::std::endl;
    }
}

int main(){
	assignment2();
	assignment5();
}
