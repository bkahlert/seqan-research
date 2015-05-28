
#include <seqan/index.h>
#include <seqan/find.h>
#include <iostream>


using namespace std;
using namespace seqan;



int main(){
	typedef StringSet<String<int>> x; 
	x y;
	String<int>z;
	appendValue(z,1);
	appendValue(z,2);
	appendValue(z,3);
	String<int>p = z;
	appendValue(y,p);
	appendValue(z,4);
	appendValue(y,z);
	Index <x> index(y);
	Finder <Index<x>> finder(y);
	Pattern <String<int>> pattern(p);
	while(find(finder,pattern)){
		cout <<"["<<beginPosition(finder)<<";"<<endPosition(finder)<<")\t"<<infix(finder)<<endl;
	}
	return 0;
}