#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <iostream>
#include <seqan/file.h>
#include <seqan/modifier.h>

using namespace seqan;

Dna getRevCompl(Dna const & nucleotide)
{
    if(nucleotide == (Dna)'A')
        return (Dna)'T';
    if(nucleotide == (Dna)'T')
        return (Dna)'A';
    if(nucleotide == (Dna)'C')
        return (Dna)'G';
    return (Dna)'C';
}
	//************************************************************
	// Function belongs to Assignment 4
// Function to print simple alignment between two sequences with the same length
// .. for two sequences of the same type
template <typename TText, typename TRead>
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
	// Assignment 1
    DnaString genome = "TATATACGCGCGAGTCGT";
    DnaString revComplGenome;
    resize(revComplGenome,length(genome));
	std::cout << length(genome) << std::endl;
	for(unsigned i=0; i<length(genome); ++i) {
		revComplGenome[i] = getRevCompl(genome[length(genome)-1-i]);
	}
	std::cout << "my result: " << std::endl;
	std::cout << revComplGenome << std::endl;
    // Your code snippet

    // And to check if your output is correct, 
    // use the given SeqAn function reverseComplement(),
    // which modifies the sequence in-place
    reverseComplement(genome);
    std::cout << genome << std::endl;
	std::cout << "***********************************************************************" << std::endl;
	std::cout << "Assignment 2" << std::endl;

	String<Dna5> nucleotides = "AGTCGTGNNANCT";
    String<Dna5String> selected;
	resize(selected,2);
	int c1 = 0;
	int c2 = 0;

    // Append all elements of nucleotides, apart of Gs, 
    // to the list `selected` 
    for (unsigned i = 0; i < length(nucleotides); ++i){
		if(nucleotides[i] < (Dna) 'G'){appendValue(selected[0], nucleotides[i]);++c1;}
		else if(nucleotides[i] > (Dna) 'G'){appendValue(selected[1], nucleotides[i]);++c2;}
        //appendValue(selected, nucleotides[i]);
    }
    std::cout << "Selected nucleotides: " << selected[0] << selected[1] << std::endl;

	std::cout << "***********************************************************************" << std::endl;
	std::cout << "Assignment 3" << std::endl;

	String<Dna5String> tmpString;
	resize(tmpString,3);
	tmpString[0] = "ATATANGCGT";
	tmpString[1] = "AAGCATGANT";
	tmpString[2] = "TGAAANTGAC";
	String<Dna5> bsp = "GATGCATGAT";
	String<Dna5String> lesser;
	String<Dna5String> greater;

	for (unsigned i = 0; i < length(tmpString); ++i){
		Lexical <> comparator(tmpString[i], bsp);
		if(isLess(comparator)){appendValue(lesser, tmpString[i]);}
		else{appendValue(greater, tmpString[i]);}
        //appendValue(selected, nucleotides[i]);
    }
	std::cout << tmpString << std::endl;

	std::cout << "***********************************************************************" << std::endl;
	std::cout << "Assignment 4" << std::endl;

	 // We have given a genome sequence
    Dna5String genome2 = "ATGGTTTCAACGTAATGCTGAACATGTCGCGT";
    // A read sequence
    Dna5String read = "TGGTNTCA";
    // And the begin position of a given alignment between the read and the genome
    unsigned beginPosition = 1;
	
	Infix<String<Dna5> >::Type genomeFragment = infix(genome2, beginPosition, (beginPosition+length(read)));
    //Dna5String genomeFragment;       
    // We have to create a copy of the corresponding fragment of the genome, where the read aligns to
    //for (unsigned i = 0; i < length(read); ++i){
    //    appendValue(genomeFragment, genome[beginPosition+i]);
    //}
    // Call of our function to print the simple alignment
    printAlign(genomeFragment, read);
  

    return 0;
}
