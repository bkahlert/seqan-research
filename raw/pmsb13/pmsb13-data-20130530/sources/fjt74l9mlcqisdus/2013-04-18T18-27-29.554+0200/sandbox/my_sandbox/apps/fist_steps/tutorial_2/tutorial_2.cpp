#include <iostream>
#include <seqan/file.h>
#include <seqan/sequence.h>
#include <seqan/score.h>

template <typename T_tex,typename T_pat>
int computeLocalScore(T_tex const & subText, T_pat const & pattern)
{
    int localScore = 0;
    for (unsigned i = 0; i < seqan::length(pattern); ++i)
        if (subText[i] == pattern[i])
            ++localScore;
    
    return localScore;
}


template <typename T_tex,typename T_pat>
seqan::String<int> computeScore(T_tex const & text, T_pat const & pattern)
{
    seqan::String<int> score;
	seqan::resize(score, seqan::length(text)-seqan::length(pattern)+1, 0);

    for (unsigned i = 0; i < seqan::length(text) - seqan::length(pattern) + 1; ++i)
        score[i] = computeLocalScore(infix(text, i, i + seqan::length(pattern)), pattern);
    
    return score;
}

void print(seqan::String<int> const & score){
	for (unsigned i = 0; i < seqan::length(score); ++i)
        std::cout << score[i] << " ";
    std::cout << std::endl;
}

template <typename T_typ>
void print(T_typ const & score){
        std::cout << score<<std::endl;
}

int main()
{
    seqan::String<char> text = "This is an awesome tutorial to get to now SeqAn!";
    seqan::String<char> pattern = "tutorial";
    seqan::String<int> score = computeScore(text, pattern);
	print(score);
	print(text);
    system("Pause");
    return 0;
}