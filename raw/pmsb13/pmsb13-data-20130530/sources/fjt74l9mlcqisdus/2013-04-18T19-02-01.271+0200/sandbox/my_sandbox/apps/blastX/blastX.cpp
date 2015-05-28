#include "own_functions.h"

struct Values{
	String<char> fasta_file;
	//Values() : fasta_file("database.fa") {}
	String<char> fastq_file;
	//Values() : fastq_file("reads.fq") {}
};

int main(int argc, char const ** argv){
	Values comVal; 
	if (argc>=2){
		if (PARSE_ARGUMENTS(argc,argv,comVal))
			return 1;
	}
	return 0;
}