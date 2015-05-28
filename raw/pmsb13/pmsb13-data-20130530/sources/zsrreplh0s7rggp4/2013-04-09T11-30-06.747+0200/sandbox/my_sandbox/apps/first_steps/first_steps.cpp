#include <iostream>
#include <seqan/file.h>
#include <seqan/sequence.h>
#include <seqan/score.h>

template <typename TText>
int computeLocalScore(TText const & subText, seqan::String<seqan::AminoAcid> const & pattern)
{
    int localScore = 0;
    for (unsigned i = 0; i < seqan::length(pattern); ++i)
        localScore += seqan::score(seqan::Blosum62(), subText[i], pattern[i]);
    
    return localScore;
}

seqan::String<int> computeScore(seqan::String<char> const & text, seqan::String<char> const & pattern)
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

void seqprint(seqan::String<int> const & score)
{
  for (unsigned i = 0; i < seqan::length(score); ++i)
        std::cout << score[i] << " ";
    std::cout << std::endl;
}


int main()
{
    seqan::String<char> text = "This is an awesome tutorial to get to now SeqAn!";
    seqan::String<char> pattern = "tutorial";
    seqan::String<int> score = computeScore(text, pattern);
    
    seqprint(score);
    
    return 0;
}