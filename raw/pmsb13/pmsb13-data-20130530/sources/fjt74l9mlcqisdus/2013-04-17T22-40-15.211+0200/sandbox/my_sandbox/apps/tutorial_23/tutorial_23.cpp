#include <seqan/sequence.h>
#include <vector>

using namespace seqan;
using namespace std;

class ami{
public:
	vector<char> Amino;
	vector<int> gruppe;
};


void lala(vector<ami> & x){
	ami eins;
	append(eins.Amino,"Y");
	append(eins.gruppe,1);
	append(x,eins);
	ami zwei;
	append(zwei.Amino,"A");
	append(zwei.gruppe,2);
	append(x,zwei);
}

int main(){
	vector<ami>x;
	lala(x);
	cout << x.size()<<endl;
	return 0;
}