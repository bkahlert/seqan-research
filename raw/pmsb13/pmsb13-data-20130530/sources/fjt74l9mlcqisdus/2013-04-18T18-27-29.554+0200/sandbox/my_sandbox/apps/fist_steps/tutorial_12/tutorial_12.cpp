#include <seqan/seq_io.h>

using namespace std;

int main(int argc,char const ** argv){

	if (argc>=2){
	
		seqan::FaiIndex faiIndex;
		int res = build(faiIndex, argv[1]);
		if (res != 0) cerr << "ERROR: Could not build the index!\n";	

		int res2 = write(faiIndex);
		if (res2 != 0) cerr << "ERROR: Could not write the index to file!\n";
	}
	else return 1;
	return 0;	
}



