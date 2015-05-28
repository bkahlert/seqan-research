#include <seqan/cluster.h>
#include <seqan/structs/PerformanceSample.h>

using namespace seqan;
using namespace std;

int clusterSearch(String<PerformanceSample> & performance, int argc, char const ** argv) {	
	in res = cluster(performance, argc, argv);
	return res;
}