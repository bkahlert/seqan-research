#include <seqan/seq_io.h>

using namespace seqan;
using namespace std;

int main(int argc, char const ** argv) {
    
    seqan::FaiIndex faiIndex;
    if (build(faiIndex, argv[1],argv[2])!= 0)
    	std::cerr << "ERROR: Could not build the index!\n";

}