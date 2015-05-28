#include <iostream>
#include <seqan/sequence.h>  // CharString, ...
#include <seqan/file.h>      // to stream a CharString into cout

#include <seqan/seq_io.h>

int main(int, char const **)
{
 	seqan::SequenceStream seqStream("C:\\blatemp\\example.fasta");
	
	if (!isGood(seqStream))
    {
        std::cout << "ERROR: Could not open the file.\n";
        return 1;
    }

	seqan::CharString id;
    seqan::Dna5String seq;
	
	while( !readRecord(id, seq, seqStream) ){
		std::cout << id << '\n' << seq << '\n\n' ;
	}

	seqan::CharString mySeqanString = "Done.";
    std::cout << mySeqanString << std::endl;
	return 1;
}