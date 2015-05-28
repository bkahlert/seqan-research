#include <iostream>
#include <seqan/sequence.h>
#include <seqan/file.h>

using namespace seqan;

// Function to print simple alignment between two sequences with the same length
// .. for two sequences of the same type
template <typename T1, typename T2>
void printAlign(T1 const & genomeFragment, T2 const & read)
{
        std::cout <<  "Alignment " << std::endl;
        std::cout << "  genome : ";
        std::cout << genomeFragment << std::endl;
        std::cout << "  read   : ";
        std::cout << read << std::endl;
}

int main()
{
    // We have given a genome sequence
    Dna5String genome = "ATGGTTTCAACGTAATGCTGAACATGTCGCGT";
    // A read sequence
    Dna5String read = "TGGTNTCA";
    // And the begin position of a given alignment between the read and the genome
    unsigned beginPosition = 1;

        
    // We have to create a copy of the corresponding fragment of the genome, where the read aligns to
	Infix<Dna5String>::Type inf = infix(genome, beginPosition, beginPosition + length(read));
    // Call of our function to print the simple alignment
    printAlign(inf, read);
  
    return 0;
}
/*#include <iostream>
#include <seqan/align.h>

using namespace seqan;

int main()
{
	typedef String<char> TSequence;
	typedef Align<TSequence,ArrayGaps> TAlign;
	typedef Row<TAlign>::Type TRow;

	TSequence seq1= "ACGTCACCTC";
	TSequence seq2= "ACGGGCCTATC";
	TAlign align;
	resize(rows(align),2);
	assignSource(row(align,0),seq1);
	assignSource(row(align,1),seq2);

	TRow &row1=row(align,0);
	TRow &row2=row(align,1);

	std::cout << align << std::endl;

	insertGap(row1,5);
	insertGaps(row1,2,2);

	insertGaps(row2,9,2);

	std::cout << align << std::endl;

	int count=0;
	typedef Iterator<TRow>::Type TRowIterator;
    TRowIterator it = begin(row1);
    TRowIterator itEnd = end(row1);
    for(; it != itEnd; ++it)
    {
        if(isGap(it))
            count++;
        
    }
    std::cout << "total count of gaps: " << count << std::endl;


}*/