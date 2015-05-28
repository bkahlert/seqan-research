#include <seqan/sequence.h>

using namespace seqan;
using namespace std;


int main(){
	StringSet<String<int>> x;
	resize(x,3);
	String<int>y;
	append(y,1);
	append(y,2);
	append(y,3);
	x[0]=y;
	String<int>z;
	append(y,4);
	append(y,5);
	append(y,6);
	x[1]=z;
	cout<<length(x)<<endl;
	cout<<x[1][0]<<endl;
	return 0;
}