#include <seqan/clusterSearch.h>
#include <seqan/structs/PerformanceSample.h>

using namespace seqan;
using namespace std;

int main(int argc, char const ** argv) {
	cout << argc << endl;
	cout << length(argv) << endl;
	String<PerformanceSample> performance;
	
	int res = clusterSearch(performance, argc, argv);
	return res;
}