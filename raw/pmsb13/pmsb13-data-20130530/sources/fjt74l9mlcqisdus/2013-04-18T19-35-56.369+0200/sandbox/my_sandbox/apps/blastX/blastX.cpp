#include "own_functions.h"
#include <iostream>

using namespace std;
using namespace seqan;

int main(int argc, char const ** argv){
 
	if (argc>=2){
		if (PARSE_ARGUMENTS(argc,argv))
			return 1;
	}
	return 0;
}
