#include <seqan/clusterSearch.h>
#include <seqan/structs/PerformanceSample.h>

using namespace seqan;
using namespace std;

int main(int argcm, char const ** argmv) {
	String<PerformanceSample> performance;
	
	int res = clusterSearch(performance, argc, argv);
	return res;
}