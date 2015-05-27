// ==========================================================================
//                                   my_app
// ==========================================================================
// Copyright (c) 2006-2011, Knut Reinert, FU Berlin
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


#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <iostream>
using namespace seqan;




template <typename TAlphabet>
void showAllLetterOfMyAlphabet(TAlphabet const &)
{
    typedef typename Size<TAlphabet>::Type TSize;
    TSize alphSize = ValueSize<TAlphabet>::VALUE;
    for (TSize i = 0; i < alphSize; ++i)
        std::cout << i << ',' << TAlphabet(i) << "  ";
    std::cout << std::endl;
}

int main(int, char **)
{
	seqan::String<AminoAcid> str = "MQDRVKRPMNAFIVWSRDQRRKMALEN";
	PeptideIterator it = begin(str);
	PeptideIterator itEnd = end(str);

	
	//absFreqOfLettersInSeq(str,it);

	while (it != itEnd)
	{
		::std::cout << *it;
		if(*it=='R')
		{
			value(it)='A';
		}
		++it;
	}
	::std::cout << ::std::endl;

	goBegin(it);
	//it = begin(str);

	while (it != itEnd)
	{
		::std::cout << *it;
		++it;
	}
	::std::cout << ::std::endl;

	

		
	//String<AminoAcid> str = "MQDRVKRPMNAFIVWSRDQRRKMALEN";

    //typedef Iterator<TAminoAcidString>::Type TIter;
    //((TIter itEnd = end(str);
    //for (TIter it = begin(str); it != itEnd; goNext(it))
    //{
      //  if (value(it) == 'R') value(it) = 'A';
        //std::cout << value(it) << ',';
    //}
    //std::cout << std::endl;

	//PeptideIterator it1 = begin(str);

	//absFreqOfLettersInSeq(it);

	//showAllLetterOfMyAlphabet(AminoAcid());
//   showAllLetterOfMyAlphabet(Dna());
//	 showAllLetterOfMyAlphabet(Dna5());
//	showAllLetterOfMyAlphabet(Rna5());
    return 0;

}