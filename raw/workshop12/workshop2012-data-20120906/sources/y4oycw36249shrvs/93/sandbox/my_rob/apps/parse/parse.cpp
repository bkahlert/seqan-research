#include <iostream>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>


int main(){

	seqan::CharString id;
    seqan::Dna5String seq;
    seqan::SequenceStream seqStream("exemple.fa");
    readRecord(id, seq, seqStream);
      if (!isGood(seqStream)){
        std::cerr << "ERROR: Could not open the file.\n";
        return 1;
      }
	while (!atEnd(seqStream)){

      if (readRecord(id, seq, seqStream) != 0){
        std::cerr << "ERROR: Could not read from example.fa!\n";
        return 1;
      }
	

      std::cout << id << '\t' << seq << '\n';
    }
    return 0;
}