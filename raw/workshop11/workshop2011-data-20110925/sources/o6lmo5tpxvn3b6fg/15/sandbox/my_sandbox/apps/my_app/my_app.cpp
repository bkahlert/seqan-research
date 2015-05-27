
#include <seqan/sequence.h>  // CharString, ...
#include <seqan/file.h>      // to stream a CharString into cout

using namespace seqan
int main(int, char **) {
	seqan::String<seqan::DNA> testseq = "GATTACA";
    std::cout << "Hello World!" << std::endl;
    seqan::CharString mySeqanString = "Hello SeqAn!";
    std::cout << mySeqanString << std::endl;
    return 1;
}