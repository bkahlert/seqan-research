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
	auto y = x[0].Amino;
	return 0;
}