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

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/misc/misc_cmdparser.h>

#include "my_app.h"

using namespace seqan;

// Program entry point
int main(int argc, char const ** argv) {
	String<AminoAcid> aaString = "MQDRVKRPMNAFIVWSRDQRRKMALEN";

	Iterator<String<AminoAcid> >::Type it = begin(aaString);

	while (atEnd(it)) {
		if (*it == 'A') {
			*it = 'R';
		}
		it++;
	}


	typedef Size<String<AminoAcid> >::Type TSize;
	String<TSize> counterString;
	TSize alphabetSize = ValueSize<AminoAcid>::VALUE;
	resize(counterString, alphabetSize);

	it = begin(aaString);

	while (atEnd(it)) {
		value(counterString, ordValue(value(it))) += 1;
		it++;
	}


	Iterator<String<TSize> >::Type it2 = begin(counterString);
	TSize pos = 0;
	while (atEnd(it2)) {
		std::cout << AminoAcid(pos++) << " " << value(it2) << std::endl;
		it2++;
	}


	std::cout << aaString << std::endl;
}