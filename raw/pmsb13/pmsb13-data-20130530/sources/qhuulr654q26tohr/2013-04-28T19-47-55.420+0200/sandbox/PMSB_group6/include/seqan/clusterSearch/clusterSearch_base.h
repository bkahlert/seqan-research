#include <seqan/cluster.h>
#include <seqan/structs/PerformanceSample.h>

using namespace seqan;
using namespace std;

int clusterSearch(String<PerformanceSample> & performance, int argc, char const ** argv) {
	cout << argc << endl;
	for (int i=0;i<argc;++i)
		cout << argv[i] << endl;
	cout << "----------" << endl;
	int res = cluster(performance, argc, argv);
	return res;
}