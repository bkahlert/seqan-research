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
	vector<ami>x;
	ami eins;
	eins.Amino = "ARNDCEQGHILKMFPSTWYV";
	x.push_back(eins);
	cout << x[0]<<endl;

	return 0;
}