#include <seqan/clusterSearch.h>
#include <seqan/structs/PerformanceSample.h>

using namespace seqan;
using namespace std;

int main(int argc, char const ** argv) {
	String<PerformanceSample> performance;
	
	int res = clusterSearch(performance, argc, argv);
	
	for (int i=0; i<length(performance);++i) {
		cout << performance[i].name << endl;
		cout << '\t' << performance[i].getTime() << endl << endl;
	}
	return res;
}