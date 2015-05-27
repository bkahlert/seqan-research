// ==========================================================================
//                                motifSearch
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
// Author: Vishalini Vimalakanthan <vishalini.v@fu-berlin.de>
// ==========================================================================

#include "motifHeader.h"

using namespace seqan;

//Read sequences from fasta-file into a StringSet

template<typename TStringSet, typename TId>
void 
readSequences(CharString &filename, TStringSet &seqs, StringSet<TId> &seqIDs)
{

	typedef String<Dna> TSequence;
	typedef CharString TId;
	typedef StringSet<TSequence> TStringSet;

	MultiSeqFile multiSeqFile;
	if (!open(multiSeqFile.concat, toCString(filename), OPEN_RDONLY))
	    {
	        std::cerr << "Failed to open " << filename << " file." << std::endl;
       
	    }

    AutoSeqFormat format;
	guessFormat(multiSeqFile.concat, format);
	split(multiSeqFile, format);

	unsigned seqCount = length(multiSeqFile);

	reserve(seqs, seqCount, Exact());
	reserve(seqIDs, seqCount, Exact());

	TSequence seq;
	TId seqid;

	for (unsigned i = 0; i < seqCount; ++i)
	{

		assignSeq(seq, multiSeqFile[i], format);   
		assignSeqId(seqid, multiSeqFile[i], format);   

		appendValue(seqs, seq, Generous());
		appendValue(seqIDs, seqid, Generous());
	}


}

int main(int argc, char const ** argv)
{


	if(argc!=2)
	{

			std::cerr << "Error: Too many arguments"<<std::endl
					  << "USAGE: " << argv[0] << " FILE" << std::endl;

		    return 1;
	}

	typedef String<Dna> TSequence;
	typedef CharString TId;
	typedef StringSet<TSequence> TStringSet;
	
	CharString filename = argv[1];
	TStringSet dataset;
	StringSet<TId> seqIDs;


	//read sequences
	readSequences(filename, dataset, seqIDs);

	// input values 
	double const epsilon = 0.16;
	const int motifLength = 13;   //actually minimum length of the motif
	typedef unsigned int TSize;
	TSize xDrop = 3; 
	
	MotifFinder<Dna, MotifNewAlgorithm> finder(epsilon, motifLength, xDrop);
	
	findMotif(finder,dataset,Oops());


	return 0;

}

