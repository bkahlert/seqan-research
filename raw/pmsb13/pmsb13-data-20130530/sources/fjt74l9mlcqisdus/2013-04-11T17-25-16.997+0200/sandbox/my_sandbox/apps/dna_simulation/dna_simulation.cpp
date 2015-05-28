#include <iostream>
#include "my_header.h"

struct Values{
	unsigned DNA_LENGTH;

};


int main(int argc, char const ** argv)
{
	// take the comandline arguments and save them in ...
	Values comVal; 
	if (argc>=2){
		if (PARSE_ARGUMENTS(argc,argv,comVal))return 1;
	}

	 
	return 0;
}