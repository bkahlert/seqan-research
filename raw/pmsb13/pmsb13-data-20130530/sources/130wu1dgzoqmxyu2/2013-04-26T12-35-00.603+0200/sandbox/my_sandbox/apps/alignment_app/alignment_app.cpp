/*#include <iostream>
#include <seqan/align.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/index.h>

using namespace seqan;
using namespace std;

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

int seqIO_assignment1(const char* file){
    seqan::CharString id;
    seqan::Dna5String seq;
    seqan::CharString qual;

    seqan::SequenceStream seqStream(file);

    if (!isGood(seqStream))
    {
        std::cerr << "ERROR: Could not open the file.\n";
        return 1;
    }

    while( !atEnd(seqStream) ){
        if (readRecord(id, seq, qual, seqStream) != 0)
        {
            std::cerr << "ERROR: Could not read from example.fa!\n";
            return 1;
        }
        std::cout << id << '\t' << seq << " " << qual << '\n';
    }
}

void assignmentIndex2(){
	String<char> text = "TTATTAAGCGTATAGCCCTATAAATATAA";
	String<char> pattern = "TATAA";
	Index<String<char>, IndexEsa<> > esaIndex(text);
	Finder<Index<String<char>, IndexEsa<> > > esaFinder(esaIndex);

	while(find(esaFinder, pattern)){
		std::cout<<position(esaFinder)<<std::endl;
	}
}

int main(int argc, char const** args){
	/*if (argc != 2){
		std::cerr << "Usage: " << args[0] << " DATEINAME" << std::endl;
		return 1;
	}
	seqIO_assignment1(args[1]);
	assignmentIndex2();

}
