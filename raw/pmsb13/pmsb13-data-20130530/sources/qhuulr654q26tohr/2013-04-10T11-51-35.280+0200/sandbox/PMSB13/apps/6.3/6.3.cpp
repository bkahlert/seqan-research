#include <iostream>
#include <seqan/align.h>
#include <seqan/graph_msa.h>

using namespace seqan;
using namespace std;

int main()
{
    typedef String<Rna> TSequence;                 // sequence type
    typedef Align<TSequence,ArrayGaps> TAlign;      // align type
    typedef Row<TAlign>::Type TRow;                 // gapped sequence type
    
    TSequence seq1 = "AAGUGACUUAUUG";
    TSequence seq2 = "AGUCGGAUCUACUG";
    
    TAlign align;
    resize(rows(align), 2);
    assignSource(row(align,0),seq1);
    assignSource(row(align,1),seq2);
    
    int score=globalAlignment(align,Score<int,MyersHirschberg());
    cout << align << endl;
    cout << score << endl;
    
    typedef Iterator<TRow, Rooted>::Type TIt1;
    typedef Iterator<Rows<TAlign>::Type, Rooted>::Type TIt2;
    TIt2 it2 = begin(rows(align));
    TIt1 it1;
    while(!atEnd(it2)){
	it1 = begin(row(align,position(it2)));
	cout << "row " << position(it2);
	while(!atEnd(it1)){
	    if (isGap(*it2,position(it1)))
		cout << position(it1) << " ";
	    goNext(it1);
	}
	cout << endl;
	goNext(it2);
    }
    
}