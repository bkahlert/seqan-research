#include <iostream>
#include <seqan/align.h>

using namespace seqan;
using namespace std;

int main()
{
    typedef String<char> TSequence;                 // sequence type
    typedef Align<TSequence,ArrayGaps> TAlign;      // align type
    typedef Row<TAlign>::Type TRow;                 // gapped sequence type
    
    TSequence seq1 = "CDFGDC";
    TSequence seq1_mod = "AACDFGDC";
    TSequence seq2 = "CDEFGAHGC";
    
    TAlign align;
    resize(rows(align), 2);
    assignSource(row(align,0),seq1);
    assignSource(row(align,1),seq2);
    
    std::cout << align;
    TRow &row1 = rows(align)[0];
    TRow &row2 = row(align,1);
    insertGap(row1,2);
    insertGaps(row1,5,2);
    std::cout << align;
    std::cout << length(row(align,0)) << std::endl;
    
    
    
    std::cout << std::endl << "Before clipping:\n" << align;
    setClippedBeginPosition(row1,1);
    setClippedEndPosition(row1,7);
    setClippedBeginPosition(row2,1);
    setClippedEndPosition(row2,7);
    std::cout << std::endl << "After clipping:\n" << align;
    
    std::cout << std::endl << "pos1: ";
    for(unsigned i = 0; i < length(row1); ++i)
	std::cout << row(align, 0)[i] << ",";
    
    insert(seq1,0,"AA");
    //assignSource(row(align,0),seq1_mod);
    
    std::cout << align;
    
    std::cout << std::endl << "pos1_mod: ";
    for(unsigned i = 0; i < length(row1); ++i)
	std::cout << row(align, 0)[i] << ",";
    
    std::cout << std::endl << "ViewToSource1: ";
    for(unsigned i = 0; i < length(row1); ++i)
	std::cout << toSourcePosition(row1, i) << ",";
    
    std::cout << std::endl << "ViewToSource2: ";
    for(unsigned i = 0; i < length(row2); ++i)
	std::cout << toSourcePosition(row2, i) << ",";
    std::cout << std::endl;
    
    std::cout << std::endl << "SourceToView1: ";
    for(unsigned i = 0; i < length(source(row1)); ++i)
	std::cout << toViewPosition(row1, i) << ",";
    
    std::cout << std::endl << "SourceToView2: ";
    for(unsigned i = 0; i < length(source(row2)); ++i)
	std::cout << toViewPosition(row2, i) << ",";
    std::cout << std::endl;
}