#include <iostream>
#include <seqan/sequence.h>
#include <seqan/file.h>

using namespace seqan;

// Function to print simple alignment between two sequences with the same length
// .. for two sequences of the same type
template <typename TText>
void printAlign(TText const & genomeFragment, TText const & read)
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

    Dna5String genomeFragment;       
    // We have to create a copy of the corresponding fragment of the genome, where the read aligns to
    for (unsigned i = 0; i < length(read); ++i){
        appendValue(genomeFragment, genome[beginPosition+i]);
    }
    // Call of our function to print the simple alignment
    printAlign(genomeFragment, read);
    String <int> test1;
    appendValue(test1, 1337);
    int test3 = 102;
    int * test2 = *test3;
    std::cout<<test2<<test1[0]<<'\n';
    return 0;
}
