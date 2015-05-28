#include "own_functions.h"

class Variable{
public:
	String<char>fasta_file;
};

int main(int argc, char const ** argv){
	Variable comVal;
	if (argc>=2){
		if (PARSE_ARGUMENTS(argc,argv,comVal))
			return 1;
	}
	return 0;
}