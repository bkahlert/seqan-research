#include <iostream>
#include "my_header.h"



int main(int argc, char const ** argv)
{
	// take the comandline arguments and save them in ...
	if (argc>=2){
		if (PARSE_ARGUMENTS(argc,argv))return 1;
	}

	 
	return 0;
}