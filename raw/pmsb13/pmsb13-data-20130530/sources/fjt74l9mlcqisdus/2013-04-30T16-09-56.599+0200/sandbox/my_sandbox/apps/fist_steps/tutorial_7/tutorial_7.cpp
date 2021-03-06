#include <iostream>
#include <seqan/sequence.h>
#include <seqan/file.h>

using namespace seqan;

// Function to print simple alignment between two sequences with the same length
// .. for two sequences of the same type
template <typename TText,typename TRead>
void printAlign(TText const & genomeFragment, TRead const & read)
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

	Infix<Dna5String>::Type genomeFragment= infix(genome,beginPosition,beginPosition+length(read));       
    // We have to create a copy of the corresponding fragment of the genome, where the read aligns to
    // Call of our function to print the simple alignment
    printAlign(genomeFragment, read);
	system("Pause");
    return 0;
}