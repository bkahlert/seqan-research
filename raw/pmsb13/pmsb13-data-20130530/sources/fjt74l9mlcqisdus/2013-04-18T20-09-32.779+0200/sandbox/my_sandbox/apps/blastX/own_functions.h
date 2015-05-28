
#ifndef SANDBOX_MY_SANDBOX_APPS_BLASTX_OWN_FUNCTIONS_
#include <iostream>
#include <seqan/arg_parse.h>

using namespace std;
using namespace seqan;

class Variable{
public:
	String<char> fasta_file ;
	Variable() : fasta_file("fasta.fasta") {}
	String<char> fastq_file ;
	Variable() : fasta_file("fastq.fastq") {}
};

int PARSE_ARGUMENTS(int argc,char const ** argv,Variable & comVal);


#define SANDBOX_MY_SANDBOX_APPS_BLASTX_OWN_FUNCTIONS_
#endif
