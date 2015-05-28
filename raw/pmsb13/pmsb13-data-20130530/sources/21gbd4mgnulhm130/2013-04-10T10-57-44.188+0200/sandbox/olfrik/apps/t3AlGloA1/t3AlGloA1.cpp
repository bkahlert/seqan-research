// ==========================================================================
//                                 t3AlGloA1
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

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/arg_parse.h>
/*
Compute a global alignment between the DNA sequences "AAATGACGGATTG" "AGTCGGATCTACTG" using the Gotoh algorithm with the following scoring parameters: match = 4, mismatch = -2, gap open = -4 and gap extend = -2. Store the alignment in an Align object and and print it together with the score.
*/


#include <iostream>
#include <seqan/align.h>

using namespace seqan;

int main()
{
    typedef String<char> TSequence;                 // sequence type
    typedef Align<TSequence,ArrayGaps> TAlign;      // align type
    TSequence seq1 = "AAATGACGGATTG";
    TSequence seq2 = "CAGTCGGATCTACTG";

    TAlign align;
    resize(rows(align), 2);
    assignSource(row(align,0),seq1);
    assignSource(row(align,1),seq2);

    int score = globalAlignment(align, Score<int,Simple>(4,-2,-4,-2));
    ::std::cout << "Score: " << score << ::std::endl;
    ::std::cout << align << ::std::endl;
    
    return 0;
}

