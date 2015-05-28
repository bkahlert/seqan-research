#include <iostream>
#include <seqan/file.h>
#include <seqan/sequence.h>

template <typename TText, typename TPattern>
int computeLocalScore(Ttext const & subText, TPattern const & pattern)
{
    int localScore = 0;
    for (unsigned i = 0; i < seqan::length(pattern); ++i)
        if (subText[i] == pattern[i])
            ++localScore;
    
    return localScore;
}


seqan::String<int> computeScore(TText const & text, TPattern const & pattern)
{
    seqan::String<int> score;
    seqan::resize(score, seqan::length(text)-seqan::length(pattern)+1, 0);

    for (unsigned i = 0; i < seqan::length(text) - seqan::length(pattern) + 1; ++i)
        score[i] = computeLocalScore(infix(text, i, i + seqan::length(pattern)), pattern);
    
    return score;
}

template <typename TScore>
void seqprint(TScore const & score)
{
  for (unsigned i = 0; i < seqan::length(score); ++i)
        std::cout << score[i] << " ";
    std::cout << std::endl;
}


int main()
{
    TText text = "This is an awesome tutorial to get to now SeqAn!";
    TPattern pattern = "tutorial";
    TScore score = computeScore(text, pattern);
    
    seqprint(score);
    
    return 0;
}