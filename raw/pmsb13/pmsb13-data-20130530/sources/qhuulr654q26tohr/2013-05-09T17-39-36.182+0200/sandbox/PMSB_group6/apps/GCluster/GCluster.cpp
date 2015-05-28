#include <iostream>
#include <seqan/sequence.h>

#include <ginger/GCluster.h>

using namespace seqan;
using namespace std;

int main(int argc, char const ** argv) {
	String<PerformanceSample> performance;
	
	int res = cluster(performance, argc, argv);
	if (res!=0) {
		cerr << "Error in function cluster" << endl;
		return 1;
	}
	
	return 0;
}