#include <iostream>
#include <seqan/sequence.h>
#include <seqan/file.h>

using namespace seqan;
// Function to print simple alignment between two sequences with the same length
template <typename TText1, typename TText2>
void printAlign(TText1 const & genomeFragment, TText2 const & read)
{
    std::cout <<  "Alignment " << std::endl;
    std::cout << "  genome : " << genomeFragment << std::endl;
    std::cout << "  read   : " << read << std::endl;
}

int main(int, char const **)
{
    // Build reads and genomes
    DnaString chr1 = "TATAATATTGCTATCGCGATATCGCTAGCTAGCTACGGATTATGCGCTCTGCGATATATCGCGCTAGATGTGCAGCTCGATCGAATGCACGTGTGTGCGATCGATTAGCGTCGATCATCGATCTATATTAGCGCGCGGTATCGGACGATCATATTAGCGGTCTAGCATTTAG";
    
    // Build List containing all reads
    typedef String<DnaString> TDnaList;
    TDnaList readList;
    resize(readList, 4);
    readList[0] = "TTGCTATCGCGATATCGCTAGCTAGCTACGGATTATGCGCTCTGCGATATATCGCGCT";
    readList[1] = "TCGATTAGCGTCGATCATCGATCTATATTAGCGCGCGGTATCGGACGATCATATTAGCGGTCTAGCATT";
    readList[2] = "AGCCTGCGTACGTTGCAGTGCGTGCGTAGACTGTTGCAAGCCGGGGGTTCATGTGCGCTGAAGCACACATGCACA";
    readList[3] = "CGTGCACTGCTGACGTCGTGGTTGTCACATCGTCGTGCGTGCGTACTGCTGCTGACA";
    
    typedef String<int> TPosList;
    TPosList alignPosList;
    resize(alignPosList, 4);
    alignPosList[0] = 7;
    alignPosList[1] = 100;
    alignPosList[2] = 172;
    alignPosListt[3] = 272;
    
    // Append a second chromosome sequence fragment to chr1
    DnaString chr2 = "AGCCTGCGTACGTTGCAGTGCGTGCGTAGACTGTTGCAAGCCGGGGGTTCATGTGCGCTGAAGCACACATGCACACGTCTCTGTGTTCCGACGTGTGTCACGTGCACTGCTGACGTCGTGGTTGTCACATCGTCGTGCGTGCGTACTGCTGCTGACACATGCTGCTG";
    append(chr1, chr2);
    
    // Print readlist
    std::cout << " \n Read list: " << std::endl;
    for(unsigned i = 0; i < length(alignPosList); ++i)
        std::cout << alignPosList[i] << std::endl;
    
    std::cout << " \n Read list: " << std::endl;
    for(unsigned i = 0; i < length(alignPosList); ++i)
        std::cout << alignPosList[i] << std::endl;
    
    // Assume we have mapped the 4 reads to chr1 (and chr2) and now have the mapping start positions (no gaps).
    // Store the start position in a String alignPosList: 7, 100, 172, 272
    
    // Optional
    // Bisulfite conversion
    // Assume chr1 is beeing bisulfate treated: Copy chr1 to a new genome bsChr1 and exchange every 'C' with a 'T'
    DnaString bsChr1;
    
    // Print alignments of the reads with chr1 (or bsChr1) sequence using the function printAlign
    // and the positions in alignPosList.
    // To do that, you have to create a copy of the fragment in chr1 (bsChr1) that is aligned to the read.
    std::cout << " \n Print alignment: " << std::endl;
    for(unsigned i = 0; i < length(readList); ++i)
    {
        // Temporary copy of begin position (beginPosition) from alignPosList
        // of a given alignment between the read and the genome
        
        // Genome fragment
        DnaString genomeFragment;
        // We have to create a copy of the corresponding fragment of the genome, where the read aligns to
        
        
        // Call of our function to print the simple alignment
        printAlign(genomeFragment, readList[i]);
    }
    return 1;
}