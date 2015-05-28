
#include <seqan/modifier.h>
#include <seqan/sequence.h>
#include <iostream>

using namespace std;
using namespace seqan;




void f(StringSet<String<Dna>> & x){
	reverseComplement(x);

}



int main(){
	StringSet<String<Dna>> x;
	String<Dna> y;
	appendValue(y,'A');
	appendValue(y,'C');
	appendValue(y,'G');
	
	appendValue(x,y);
	
	String<Dna>z;
	appendValue(z,'T');
	appendValue(z,'A');
	appendValue(z,'A');
	
	appendValue(x,z);
	
	f(x);
	cout<<x[0]<<endl;
	cout<<x[1]<<endl;
	return 0;
}