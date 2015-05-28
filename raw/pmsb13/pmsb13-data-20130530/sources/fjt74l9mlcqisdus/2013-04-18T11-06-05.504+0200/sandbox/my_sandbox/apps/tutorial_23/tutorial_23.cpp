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
	append(z,1);
	append(z,2);
	append(z,3);
	append(z,4);
	append(x,z);
	cout<<x[0]<<endl;
	return 0;
}