#include <iostream>
#include <seqan/sequence.h>

#include <seqan/GCluster.h>

using namespace seqan;
using namespace std;

int main(int argc, char const ** argv) {
	String<PerformanceSample> performance;
	
	int res = GCluster(performance, argc, argv);
	if (res!=0) {
		cerr << "Error in function GCluster" << endl;
		return 1;
	}
	
	return 0;
}