#include <iostream>
#include <seqan/align.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/index.h>

using namespace seqan;
using namespace std;

void assignmentIndex2(){
	String<char> text = "TTATTAAGCGTATAGCCCTATAAATATAA";
	String<char> pattern = "TATAA";
	Index<String<char>, IndexEsa<> > esaIndex(text);
	Finder<Index<String<char>, IndexEsa<> > > esaFinder(esaIndex);

	while(find(esaFinder, pattern)){
		std::cout<<position(esaFinder)<<std::endl;
	}
}

int main(){

	return 0;
}
