#include <seqan/sequence.h>
#include <vector>

using namespace seqan;
using namespace std;

class ami{
public:
	String<AminoAcid> Amino;
	String<int> gruppe;
};



int main(){
	String<ami>x;
	ami eins;
	eins.Amino = "ARNDCEQGHILKMFPSTWYV";
	append(eins.gruppe,1);
	cout<<eins.gruppe<<endl;
	cout<<eins.Amino<<endl;
	append(x,eins);
	//cout<<x[0].gruppe<<endl;
	return 0;
}