#include <iostream>
#include <seqan/sequence.h>  // CharString, ...
#include <seqan/file.h>      // to stream a CharString into cout

#include <seqan/seq_io.h>

int main(int, char const **)
{
 	seqan::SequenceStream seqStream("C:\\temp\\example.fasta");
	
	seqan::CharString id;
    seqan::Dna5String seq;
	
	
	seqan::CharString mySeqanString = "Done.";
    std::cout << mySeqanString << std::endl;
	return 1;
}