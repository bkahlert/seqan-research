// Copy the code into a demo program and have a look at the result.

#include <iostream>
#include <seqan/file.h>
#include <seqan/sequence.h>

int main()
{
	// Initialization
	seqan::String<char> text = "This is an awesome tutorial to get to know SeqAn!";
	seqan::String<char> pattern = "tutorial";
	seqan::String<int> score;
	resize(score, seqan::length(text) - seqan::length(pattern) + 1);
	// Computation of the similarities
	// Iteration over the text (outer loop)
	for (unsigned i = 0; i < seqan::length(text) - seqan::length(pattern) + 1; ++i)
	{
		int localScore = 0;
		// Iteration over the pattern for character comparison
		for (unsigned j = 0; j < seqan::length(pattern); ++j)
		{
			if (text[i + j] == pattern[j])
				++localScore;
		}
		score[i] = localScore;
	}
	// Printing the result
	for (unsigned i = 0; i < seqan::length(score); ++i)
		std::cout << score[i] << " ";
	std::cout << std::endl;
	// > 1 0 1 0 0 1 0 0 0 0 0 0 1 0 0 0 0 1 0 8 0 1 0 0 0 0 2 0 1 0 0 1 0 3 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0

	return 0;
}
