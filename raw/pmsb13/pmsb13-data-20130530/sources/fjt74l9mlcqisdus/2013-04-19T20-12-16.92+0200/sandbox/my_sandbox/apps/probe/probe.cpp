
#include <seqan/modifier.h>
#include <seqan/sequence.h>
#include <iostream>

using namespace std;
using namespace seqan;




void f(StringSet<String<int>> & x){
	reverseComplement(x);

}



int main(){
	StringSet<String<int>> x;
	String<int>y;
	appendValue(y,1);
	appendValue(y,2);
	appendValue(y,3);
	
	appendValue(x,y);
	
	appendValue(y,4);
	appendValue(y,5);
	appendValue(y,6);
	
	appendValue(x,y);
	
	f(x);
	return 0;
}