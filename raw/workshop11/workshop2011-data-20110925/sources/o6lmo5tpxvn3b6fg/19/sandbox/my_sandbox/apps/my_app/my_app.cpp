#include <iostream>
#include <seqan/sequence.h>  // CharString, ...
#include <seqan/file.h>      // to stream a CharString into cout


int main(int, char **) {
	seqan::String<seqan::Dna> testseq = "GATTACA";
    std::cout << "Hello World!" << std::endl;
    std::cout << testseq << std::endl;

	seq::resize(testseq, 20);
	std::cout << testseq << std::endl;
    return 1;
}