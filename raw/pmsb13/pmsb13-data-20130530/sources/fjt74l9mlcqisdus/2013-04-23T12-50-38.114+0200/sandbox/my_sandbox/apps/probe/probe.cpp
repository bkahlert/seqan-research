
#include <seqan/index.h>
#include <seqan/find.h>
#include <iostream>


using namespace std;
using namespace seqan;

StringSet<String<int>> f(){
	StringSet<String<int>>z;
	String<int>x;
	appendValue(x,1);
	appendValue(x,2);
	appendValue(z,x);
	appendValue(z,x);

	return z;
}



int main(){
	StringSet<String<int>> k = f();
	return 0;
}