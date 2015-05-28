#include <iostream>
#include <seqan/sequence.h>  // CharString, ...
#include <seqan/file.h>      // to stream a CharString into cout

int main(int, char const **)
{
    std::cout << "Hello World!" << std::endl;
    seqan::CharString mySeqanString = "Hello SeqAn!";
    std::cout << mySeqanString << std::endl;
	return 1;}
