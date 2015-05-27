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
	cout << AA << endl;
	// get size
	cout << length(AA) << endl;
    return 0;
}