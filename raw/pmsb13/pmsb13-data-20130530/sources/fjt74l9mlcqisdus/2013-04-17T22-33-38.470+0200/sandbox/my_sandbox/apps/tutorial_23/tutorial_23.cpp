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
	append(x,eins);
	
	return 0;
}