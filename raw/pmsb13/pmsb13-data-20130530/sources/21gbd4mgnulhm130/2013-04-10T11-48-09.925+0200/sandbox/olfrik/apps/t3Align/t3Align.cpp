// ==========================================================================
//                                  t3Align
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
#include <iostream>
#include <seqan/align.h>

using namespace seqan;

/*
   Compute a semi-global alignment between the sequences 
   AAATGACGGATTG and TGGGA using the Levenshtein distance 
   and using an AlignmentGraph to store the alignment. 
   Print the score and the resulting alignment to the standard output.
 *
 */
int asmnt2()
{
    typedef String<char> TSequence;                 // sequence type
    typedef StringSet<TSequence> TStringSet;       // container for strings
    typedef StringSet<TSequence, Dependent<> > TDepStringSet;   // dependent string set
    typedef Graph<Alignment<TDepStringSet> > TAlignGraph;   // alignment graph

    TSequence seq1 = "AAATGACGGATTG";
    TSequence seq2 = "TGGGA";

    TStringSet sequences;
    appendValue(sequences,seq1);
    appendValue(sequences,seq2);

    TAlignGraph alignG(sequences);

    int score = globalAlignment(alignG, Score<int,Simple>(1,-1,-1), AlignConfig<true,false, false, true>());
    ::std::cout << "Score: " << score << ::std::endl;
    ::std::cout << alignG << ::std::endl;

    return 0;


}


/*
 *
Write a program that computes a fast global alignment between the  
Rna sequences AAGUGACUUAUUG and AGUCGGAUCUACUG using the Align data 
structure and the Levenshtein distance. Print the score and the alignment.
Additionally, output for each row of the Align object the view positions of the gaps.
 */ 
int asmnt3(){
    typedef String<Rna> TSequence;                 // sequence type
    typedef Align<TSequence, ArrayGaps> TAlign;     // align type
    typedef Row<TAlign>::Type TRow;                 // gapped sequence type
    
    TSequence seq1 = "AAGUGACUUAUUG";
    TSequence seq2 = "AGUCGGAUCUACUG";

    TAlign align;
    resize(rows(align),2);
    assignSource(row(align,0),seq1);
    assignSource(row(align,1),seq2);
     int score = globalAlignment(align,Score<int,Simple>(1,-1,-1),Hirschberg());
    ::std::cout << "Score: " << score << ::std::endl;
    ::std::cout << align << ::std::endl;

    typedef Iterator<TRow>::Type TRowIterator;
    TRowIterator it = begin(row(align,0));
    TRowIterator itEnd = end(row(align,0));
    

    
    for(; it != itEnd; ++it)
    {
        if(isGap(it)){
            std::cout << " Gap1 " << std::endl;
        }
    }
    
    return 0;
}

int main(){
   return asmnt3();
}

