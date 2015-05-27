#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/misc/misc_cmdparser.h>

#include "my_app.h"

using namespace std;
using namespace seqan;

// Program entry point
int main(int argc, char const ** argv)
{
	// create AA seq and iterate over
	typedef String<AminoAcid> TAminoAcidString;
	TAminoAcidString AA = "MQDRVKRPMNAFIVWSRDQRRKMALEN";
	cout << AA << endl;
	// get size
	int length = length(AA);
	// create iterator
	typedef Iterator<String<char> >::Type TIter;
	TIter itrStart = begin(AA);
	TIter itrEnd = end(AA);
	for(int i = itrStart; i<length; goNext(i)){
		cout << i << ", " << AA(i) << " ";
	}
    return 0;
}