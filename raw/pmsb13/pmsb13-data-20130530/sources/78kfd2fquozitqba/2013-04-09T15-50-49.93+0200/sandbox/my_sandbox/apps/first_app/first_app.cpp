#include <iostream>
#include <seqan/sequence.h>
#include <seqan/file.h>

using namespace seqan;

// Function to print simple alignment between two sequences with the same length
// .. for two sequences of the same type
template<typename TText1, typename TText2>
void printAlign(TText1 const & genomeFragment, TText2 const & read)
{
    std::cout << "Alignment " << std::endl;
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

    Infix<Dna5String>::Type genomeFragment = infix(
            genome, beginPosition, beginPosition + length(read));
    // Call of our function to print the simple alignment
    printAlign(genomeFragment, read);

    return 0;
}
