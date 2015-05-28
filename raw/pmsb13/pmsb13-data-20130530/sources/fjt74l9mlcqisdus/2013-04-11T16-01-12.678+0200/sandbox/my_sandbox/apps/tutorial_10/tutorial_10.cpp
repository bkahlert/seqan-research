#include <iostream>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

using namespace std;

int main(int argc,const char ** argv)
{
	if (argc>=2){
		seqan::CharString id;
		seqan::Dna5String seq;
		seqan::CharString qual;

		seqan::SequenceStream seqStream(argv[1]);
		if (!isGood(seqStream)){
				cerr << "ERROR: Could not open the file.\n";
				return 1;
		}
		
		
		while(!atEnd(seqStream)){
			if (readRecord(id, seq, qual,seqStream) != 0){
				cerr << "ERROR: Could not read from fq!\n";
				return 1;
			}
			cout << id << '\t' << seq <<"\t"<< qual<<endl;
		}
	}
	else return 1;
	return 0;
}