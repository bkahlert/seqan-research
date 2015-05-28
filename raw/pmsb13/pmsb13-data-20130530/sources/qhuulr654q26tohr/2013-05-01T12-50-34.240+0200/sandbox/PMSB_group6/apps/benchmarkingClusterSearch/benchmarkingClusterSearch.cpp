#include <seqan/clusterSearch.h>
#include <seqan/structs/PerformanceSample.h>
#include <seqan/structs/MemorySample.h>

using namespace seqan;
using namespace std;

int main(int argc, char const ** argv) {
	String<PerformanceSample> performance;
	
	
	MemorySample ms("gesamt");
	ms.printToStdout();
	
	int res = clusterSearch(performance, argc, argv);
	ms.computeMemory();
	
	for (int i=0; i<length(performance);++i) {
		performance[i].printToStdout();
	}
	ms.printToStdout();
	return res;
}