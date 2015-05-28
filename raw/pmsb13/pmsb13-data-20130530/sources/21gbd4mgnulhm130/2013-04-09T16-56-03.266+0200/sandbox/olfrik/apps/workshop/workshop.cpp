// ==========================================================================
//                                  workshop
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
#include <iostream>
#include <seqan/file.h>
#include <seqan/modifier.h>
#include <seqan/arg_parse.h>


using namespace seqan;


// Function to print simple alignment between two sequences with the same length
template <typename TText1, typename TText2>
void printAlign(TText1 const & genomeFragment, TText2 const & read)
{
        std::cout <<  "Alignment " << std::endl;
        std::cout << "  genome : " << genomeFragment << std::endl;
        std::cout << "  read   : " << read << std::endl;
}

int main(int, char const **)
{
    // Build reads and genomes
    DnaString chr1 = "TATAATATTGCTATCGCGATATCGCTAGCTAGCTACGGATTATGCGCTCTGCGATATATCGCGCTAGATGTGCAGCTCGATCGAATGCACGTGTGTGCGATCGATTAGCGTCGATCATCGATCTATATTAGCGCGCGGTATCGGACGATCATATTAGCGGTCTAGCATTTAG";

    // Build List containing all reads
    typedef String<Dna5String> DnaList;
    DnaList readList;
    resize(readList, 4);
    readList[0] = "TTGCTATCGCGATATCGCTAGCTAGCTACGGATTATGCGCTCTGCGATATATCGCGCT";
    readList[1] = "TCGATTAGCGTCGATCATCGATCTATATTAGCGCGCGGTATCGGACGATCATATTAGCGGTCTAGCATT";
    readList[2] = "AGCCTGCGTACGTTGCAGTGCGTGCGTAGACTGTTGCAAGCCGGGGGTTCATGTGCGCTGAAGCACACATGCACA";
    readList[3] = "CGTGCACTGCTGACGTCGTGGTTGTCACATCGTCGTGCGTGCGTACTGCTGCTGACA";
    
    // Append a second chromosome sequence fragment to chr1
    DnaString chr2 = "AGCCTGCGTACGTTGCAGTGCGTGCGTAGACTGTTGCAAGCCGGGGGTTCATGTGCGCTGAAGCACACATGCACACGTCTCTGTGTTCCGACGTGTGTCACGTGCACTGCTGACGTCGTGGTTGTCACATCGTCGTGCGTGCGTACTGCTGCTGACACATGCTGCTG";
    append(chr1, chr2);
   
    // Print readlist
    std::cout << " \n Read list: " << std::endl;
    for(unsigned i = 0; i < length(readList); ++i)
        std::cout << readList[i] << std::endl;

    // Assume we have mapped the 4 reads to chr1 (and chr2) and now have the mapping start positions (no gaps).
    // Store the start position in a String: 7, 100, 172, 272
    String<unsigned> alignPosList;
    resize(alignPosList, 4);
    alignPosList[0] = 7;
    alignPosList[1] = 100;
    alignPosList[2] = 172;
    alignPosList[3] = 272;

    // Print alignments using Segment
    std::cout << " \n Print alignment using Segment: " << std::endl;
    
    Iterator<Dna5String>::Type it = begin(readList);
    for(; !atEnd(it, readList); goNext(it))
    {
    
        // Begin and end position of a given alignment between the read and the genome
        unsigned beginPosition = getValue(alignPosList);
        unsigned endPosition = beginPosition + length( getValue(readList));
        // Build infix
        Infix<DnaString>::Type genomeFragment = infix(chr1, beginPosition, endPosition);
        // Call of our function to print the simple alignment
        printAlign(genomeFragment, getValue(readList);
    }

    // Iterators :)
    // Print alignments using Iterators: Do the same as above, but use Iterators to iterate over the read list.
    // First, use Standard Iterators: Build two iterators {{{it}}} and {{{itEnd}}} to traverse {{{readList}}}.
    
    std::cout << " \n Print alignment using Standard Iterators: " << std::endl;

 
    return 1;
}
