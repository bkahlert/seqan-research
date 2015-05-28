#include "own_functions.h"






int main(int argc, char const ** argv){
	Values comVal; 
	if (argc>=2){
		if (PARSE_ARGUMENTS(argc,argv,comVal))
			return 1;
	}
	return 0;
}