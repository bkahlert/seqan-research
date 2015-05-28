// ==========================================================================
//                                 first_app
// ==========================================================================
// Copyright (c) 2006-2012, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Your Name <your.email@example.net>
// ==========================================================================

#include <iostream>
#include <seqan/file.h>
#include <seqan/sequence.h>
#include <seqan/score.h>


// ==========================================================================
// Functions
// ==========================================================================



template <typename TText, typename TPattern>
int computeLocalScore(const TText & subText, 
		      const TPattern & pattern)
{
    int localScore = 0;
    for (unsigned i = 0; i < seqan::length(pattern); ++i)
        if (subText[i] == pattern[i])
            ++localScore;
    return localScore;
}

template <typename TText>
int computeLocalScore(TText const & subText, 
		      seqan::String<seqan::AminoAcid> const & pattern)
{
    int localScore = 0;
    for (unsigned i = 0; i < seqan::length(pattern); ++i)
        localScore += seqan::score(seqan::Blosum62(), subText[i], pattern[i]);
    
    return localScore;
}





template <typename TText, typename TPattern>
seqan::String<int> computeScore(const TText & text, 
				const TPattern & pattern)
{
    seqan::String<int> score;
    seqan::resize(score, seqan::length(text), 0);
    for (unsigned i = 0; i < seqan::length(text) - seqan::length(pattern) + 1; ++i)
        score[i] = computeLocalScore(infix(text, i, i + seqan::length(pattern)), pattern);
    return score;
}

template <typename TText>
void print(const TText & text)
{
  std::cout << test << std::endl;
}

void print(const seqan::String<int> & text)
{
    for (unsigned i = 0; i < seqan::length(text); ++i)
        std::cout << text[i] << " ";
    std::cout << std::endl;
}




template <typename TText, typename TSpec>
void print(TText const & text, TSpec const & /*tag*/)
{
    print(text);
}



struct MaxOnly {};
template <typename TText>
void print(TText const & score, MaxOnly const & /*tag*/)
{
    int maxScore = score[0];
    seqan::String<int> output;
    appendValue(output, 0);
    for (unsigned i = 1; i < seqan::length(score); ++i)
    {
        if (score[i] > maxScore)
        {
            maxScore = score[i];
            clear(output);
            resize(output, 1, i);
        }
        else if (score[i] == maxScore)
            appendValue(output, i);
    }
    print(output);
}



struct GreaterZero {};
template <typename TText>
void print(TText const & score, GreaterZero const & /*tag*/)
{
    seqan::String<seqan::Pair<int> > output;
    for (unsigned i = 1; i < seqan::length(score); ++i)
        if (score[i] > 0)
            appendValue(output, seqan::Pair<int>(i, score[i]));
    for (unsigned i = 0; i < seqan::length(output); ++i)
        std::cout << "(" << output[i].i1 << "; " << output[i].i2 << ") ";
    std::cout << std::endl;
}






// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

// Program entry point.


int main()
{
    seqan::String<char> text = "This is an awesome tutorial to get to know SeqAn!";
    seqan::String<char> pattern = "tutorial";
    seqan::String<int> score = computeScore(text, pattern);
    print(text);
    print(pattern);
    print(score);
    print(score, MaxOnly());
    print(score, GreaterZero());

    // And now for a protein pattern
    seqan::String<seqan::AminoAcid> protein = "tutorial";
    seqan::String<int> proteinScore = computeScore(text, protein);
    print(text);
    print(protein);
    print(proteinScore);
    print(proteinScore, MaxOnly());
    print(proteinScore, GreaterZero());


    return 0;
}