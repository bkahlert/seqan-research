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
	append(y,4);
	append(x,y);
	//cout<<x<<endl;
	return 0;
}