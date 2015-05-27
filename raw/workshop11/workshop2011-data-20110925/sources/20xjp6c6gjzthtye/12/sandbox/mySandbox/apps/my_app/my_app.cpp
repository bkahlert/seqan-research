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
	String<char> AA = "MQDRVKRPMNAFIVWSRDQRRKMALEN";
	// create iterator
	Iterator<String<char>, Standard > itrStart = begin(AA);
	Iterator<String<char>, Standard > itrEnd = end(AA);
	// get length
	for(int i = itrStart; i < ValueSize<AA>; ++i){
		cout << i << ", " << AA(i) << " ";
	}

    return 0;
}