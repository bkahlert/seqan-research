#include <seqan/clusterSearch.h>
#include <seqan/structs/PerformanceSample.h>

using namespace seqan;
using namespace std;

int main(int argc, char const ** argv) {
	String<PerformanceSample> performance;
	
	in res = clusterSearch(performance, argc, argv);
	return res;
}