#include <iostream>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/stream.h>

using namespace seqan;
using namespace std;

int main(int argc, char const ** argv) {
    
    if (argc!=6){
	std::cerr << "USAGE: 8.2 FAIFILE FASTFILE ID START END\n";
	return 1;
    }
    
    seqan::FaiIndex faiIndex;
    if (read(faiIndex, argv[2],argv[1])!= 0){
	std::cerr << "ERROR: Could not read the index!\n";
	return 1;
    }
    
    unsigned idx = 0;
    if (!getIdByName(faiIndex, argv[3], idx))
	std::cerr << "ERROR: FAI index has no entry for "<<argv[3]<<".\n";
    
    seqan::CharString seqChr1Infix;
    unsigned start;
    lexicalCast2(start, argv[4]);
    unsigned end;
    lexicalCast2(end, argv[5]);
    if (readRegion(seqChr1Infix, faiIndex, idx, start, end) != 0)
	std::cerr << "ERROR: Could not read region.\n";
    
    cout << seqChr1Infix << endl;
    
    return 0;
}