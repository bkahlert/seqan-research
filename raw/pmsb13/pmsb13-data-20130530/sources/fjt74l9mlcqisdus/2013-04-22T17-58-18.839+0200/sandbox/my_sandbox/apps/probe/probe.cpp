
#include <seqan/index.h>
#include <seqan/find.h>
#include <iostream>


using namespace std;
using namespace seqan;



int main(){
	typedef StringSet<String<char>> x; 
	x y;
	appendValue(y,"SCHLA");
	appendValue(y,"MPELA");
	Index <x> index(y);
	Finder <Index<x>> finder(y);
	Pattern <String<char>> pattern("LA");
	while(find(finder,pattern)){
		cout <<"["<<beginPosition(finder)<<";"<<endPosition(finder)<<")\t"<<infix(finder)<<endl;
	}
	return 0;
}