#include <seqan/basic.h>
#include <seqan/file.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

using namespace seqan;

int main (int argc, char *argv[]){

	for (int i = 0; i < argc; i++ ) {
      ::std::cout<<argv[i];
   }



	return 0 ;

}