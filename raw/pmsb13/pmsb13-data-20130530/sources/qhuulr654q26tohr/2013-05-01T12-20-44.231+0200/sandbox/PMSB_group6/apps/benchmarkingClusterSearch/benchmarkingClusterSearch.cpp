#include <seqan/clusterSearch.h>
#include <seqan/structs/PerformanceSample.h>
#include <seqan/structs/MemorySample.h>

using namespace seqan;
using namespace std;

int main(int argc, char const ** argv) {
	String<PerformanceSample> performance;
	
	
	MemorySample ms("gesamt");
	cout << ms.name << endl;
	cout << '\t' << ms.virtMemory << endl;
	cout << '\t' << ms.physMemory << endl << endl;
	
	int res = clusterSearch(performance, argc, argv);
	ms.computeMemory();
	
	for (int i=0; i<length(performance);++i) {
		cout << performance[i].name << endl;
		cout << '\t' << performance[i].getTime() << endl << endl;
	}
	cout << ms.name << endl;
	cout << '\t' << ms.virtMemory << endl;
	cout << '\t' << ms.physMemory << endl << endl;
	return res;
}