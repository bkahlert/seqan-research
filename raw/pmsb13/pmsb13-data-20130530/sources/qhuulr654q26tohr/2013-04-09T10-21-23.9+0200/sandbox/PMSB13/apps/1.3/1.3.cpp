#include <iostream>
#include <seqan/file.h>
#include <seqan/sequence.h>

template <typename TText, typename TPattern>
int computeLocalScore(TText const & subText, TPattern const & pattern)
{
    int localScore = 0;
    for (unsigned i = 0; i < seqan::length(pattern); ++i)
	if (subText[i] == pattern[i])
	    ++localScore;
	
	return localScore;
}

template <typename TText, typename TPattern>
seqan::String<int> computeScore(TText const & text, TPattern const & pattern)
{
    seqan::String<int> score;
    seqan::resize(score, seqan::length(text) - seqan::length(pattern) + 1);
    
    for (unsigned i = 0; i < length(score); ++i)
	score[i] = computeLocalScore(infix(text, i, i + seqan::length(pattern)), pattern);
    
    return score;
}

template <typename TScore>
void print(TScore const & score){
    for (unsigned i = 0; i < seqan::length(score); ++i)
	std::cout << score[i] << " ";
    std::cout << std::endl;
}


int main()
{
    seqan::String<char> text = "This is an awesome tutorial to get to now SeqAn!";
    seqan::String<char> pattern = "tutorial";
    seqan::String<int> score = computeScore(text, pattern);
    print(score);
    return 0;
}
