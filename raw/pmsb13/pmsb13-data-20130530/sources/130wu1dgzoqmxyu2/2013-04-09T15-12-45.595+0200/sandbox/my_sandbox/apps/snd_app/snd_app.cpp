#include <iostream>
#include <seqan/sequence.h>
#include <seqan/file.h>

using namespace seqan;

// Function to print simple alignment between two sequences with the same length
// .. for two sequences of the same type
template <typename TFragment, typename TText>
void printAlign(TFragment const & genomeFragment, TText const & read)
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

    // Call of our function to print the simple alignment
    printAlign(infix(genome, beginPosition, length(read) + 1), read);

    return 0;
}
