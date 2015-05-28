#include <iostream>
#include <seqan/file.h>
#include <seqan/sequence.h>

template <typename TAlphabet>
int computeLocalScore(seqan::String<TAlphabet> const & subText, seqan::String<TAlphabet> const & pattern)
{
    int localScore = 0;
    for (unsigned i = 0; i < seqan::length(pattern); ++i)
	if (subText[i] == pattern[i])
	    ++localScore;
	
	return localScore;
}

template <typename TAlphabet>
seqan::String<int> computeScore(seqan::String<TAlphabet> const & text, seqan::String<TAlphabet> const & pattern)
{
    seqan::String<int> score;
    seqan::resize(score, seqan::length(text) - seqan::length(pattern) + 1);
    
    for (unsigned i = 0; i < length(score); ++i) {
	seqan::String<TAlphabet> t = infix(text, i, i + seqan::length(pattern));
	score[i] = computeLocalScore(t, pattern);
    }
    
    return score;
}

template <typename TScore>
void print(seqan::String<TScore> const & score){
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
