// ==========================================================================
//                                  t4Index
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
#include <seqan/sequence.h>
#include <seqan/index.h>

using namespace seqan;

int ass1()
{
    typedef String<Dna5> DnaString;
    typedef StringSet<DnaString> DnaStringSet;
    DnaString text = "TTATTAAGCGTATAGCCCTATAAATATAA";
    DnaStringSet patternSet;
    appendValue(patternSet, "TATA");
    appendValue(patternSet, "GC");
    Index<DnaString > index(text);
    
    Finder<Index<DnaString> > finder(index);

    //std::cout << "Searching for " << length(pattern) << " patterns in " << text << std::endl;  

    Iterator<DnaStringSet>::Type it = begin(patternSet);
    for(;!atEnd(it);++it){
       clear(finder);
        while(find(finder, getValue(it))){
            std::cout << "found " << getValue(it) << "  in pos: " << position(finder) << std::endl;
        };
    }
    return 0;
}

int ass2(){
    typedef StringSet<String<char > > CharStringSet;  
    typedef Index<StringSet<String<char> > > TIndex;
    CharStringSet text;
    appendValue(text, "How many");
    appendValue(text, " wood would");
    appendValue(text, " a woodchuck chuck?");
    TIndex index(text);
    Iterator< TIndex, TopDown< > >::Type it(index);
    CharString pattern = "wood";
    while (repLength(it) < length(pattern))
    {
        if (!goDown(it, pattern[repLength(it)])) return 0;
        unsigned endPos = std::min(repLength(it), length(pattern));
        
        std::cout << representative(it) << std::endl;
        
        if (infix(representative(it), parentRepLength(it) + 1, endPos) !=
            infix(pattern, parentRepLength(it) + 1, endPos)) {
            return 0;
        }
    }
    for (unsigned i = 0; i < length(getOccurrences(it)); ++i)
        std::cout << getOccurrences(it)[i] << std::endl;
    return 0;
}

int assq(){

    typedef Index<DnaString, IndexQGram< OneGappedShape< > > > TIndex;
    TIndex index("CATGATTACAT");
/*
    hash(indexShape(index), "AT-A");
    for (unsigned i = 0; i < length(getOccurrences(index, indexShape(index))); ++i)
        std::cout << getOccurrences(index, indexShape(index))[i] << std::endl;
*/
    
    return 0;
}

int main(){
    return assq();
}

