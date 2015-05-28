#include <seqan/sequence.h>

using namespace seqan;
using namespace std;


int main(){
	String<String<int>> x;
	String<int>y;
	append(y,1);
	append(y,2);
	append(y,3);
	append(x,y);
	String<int>z;
	append(z,4);
	append(z,5);
	append(z,6);
	append(z,7);
	append(x,z);
	cout<<length(x)<<endl;
	cout<<x[65][0]<<endl;
	return 0;
}