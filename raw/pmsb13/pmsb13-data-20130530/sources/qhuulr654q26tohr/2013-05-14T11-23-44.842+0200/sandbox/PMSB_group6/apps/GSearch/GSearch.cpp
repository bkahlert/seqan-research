#include <iostream>
#include <seqan/sequence.h>

#include <seqan/GSearch.h>

using namespace seqan;
using namespace std;

int main(int argc, char const ** argv) {
	String<PerformanceSample> performance;
	String<GMatch<int> > matches;
	
	int res = GSearch(matches, performance, argc, argv);
	if (res!=0) {
		cerr << "Error in function GSearch" << endl;
		return 1;
	}
	
	return 0;
}