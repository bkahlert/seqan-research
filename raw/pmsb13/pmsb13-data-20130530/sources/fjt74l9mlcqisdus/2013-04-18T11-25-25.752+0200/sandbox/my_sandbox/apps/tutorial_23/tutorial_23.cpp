#include <seqan/sequence.h>

using namespace seqan;
using namespace std;


int main(){
	StringSet<String<char>> x;
	x[0]="ahjdslj";
	x[1]="kewhoaf";
	cout<<x[0]<<endl;
	return 0;
}