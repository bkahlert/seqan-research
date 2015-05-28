#include <iostream>
#include <seqan/sequence.h>  // CharString, ...
#include <seqan/file.h>      // to stream a CharString into cout
#include <stdlib.h>

using namespace std;

int main(int, char const **)
{
    cout << "Hello World!" << endl;
    seqan::CharString mySeqanString = "Hello SeqAn!";
    cout << mySeqanString << endl;
	system("Pause");
	return 1;
}