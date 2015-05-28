
#include <seqan/index.h>
#include <seqan/find.h>
#include <iostream>


using namespace std;
using namespace seqan;

String<int> f(){
	String<int>x;
	appendValue(x,1);
	appendValue(x,2);
	appendValue(x,3);
	return x;
}



int main(){
	cout << f()<<endl;
	return 0;
}