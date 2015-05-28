#include <seqan/sequence.h>
#include <vector>

using namespace seqan;
using namespace std;

class ami{
public:
	vector<char> Amino;
	vector<int> gruppe;
};



int main(){
	vector<ami>x;
	ami eins;
	append(eins.Amino,"Y");
	append(eins.gruppe,1);
	append(x,eins);
	//cout<<x[0].gruppe<<endl;
	return 0;
}