
#include <seqan/index.h>
#include <seqan/find.h>
#include <iostream>


using namespace std;
using namespace seqan;





int main(){
	String<AminoAcid> x;
	appendValue(x,1);
	cout<<x<<endl; 
	return 0;
}