#include "own_functions.h"



int main(int argc, char const ** argv){
	Variable comVal;
	dafault_values(comVal);
	if (argc>=2){
		if (PARSE_ARGUMENTS(argc,argv,comVal))
			return 1;
	}
	cout<<comVal.fasta_file<<"\t"<<comVal.fastq_file<<"\t"<<comVal.seed<<
		"\t"<<comVal.size_alp<<"\t"<<comVal.numb_alp<<endl;
	return 0;
}