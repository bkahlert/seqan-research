#include <iostream>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

using namespace seqan;
using namespace std;

int main(int argc, char const ** argv) {
    
    if (argc!=3){
	std::cerr << "USAGE: FASTAFILE FAIPATH 7.6 \n";
	return 1;
    }
    
    seqan::FaiIndex faiIndex;
    if (build(faiIndex, argv[1],argv[2])!= 0){
    	std::cerr << "ERROR: Could not build the index!\n";
	return 1;
    }
    
    if (write(faiIndex)!=0){
	std::cerr << "ERROR: Could not write the index!\n";
	return 1;
    }
	
    return 0;
}