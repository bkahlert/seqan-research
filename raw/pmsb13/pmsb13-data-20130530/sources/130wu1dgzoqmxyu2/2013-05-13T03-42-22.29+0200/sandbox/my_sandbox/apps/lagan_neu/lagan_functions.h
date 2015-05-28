// ==========================================================================
//                             lagan_functions.h
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

#ifndef SANDBOX_MY_SANDBOX_APPS_LAGAN_NEU_LAGAN_FUNCTIONS_H_
#define SANDBOX_MY_SANDBOX_APPS_LAGAN_NEU_LAGAN_FUNCTIONS_H_


// ============================================================================
// Forwards
// ============================================================================

#include <seqan/seq_io.h>
#include <seqan/score.h>
#include <seqan/seeds.h>
#include <seqan/align.h>
//#include "stellar.h"

using namespace seqan;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================



bool get_global_seed_chain(seqan::String<seqan::Seed<seqan::Simple> > &seedChain,
							seqan::DnaString &seq1, seqan::DnaString &seq2,
							std::pair<unsigned, unsigned> startPos,
							std::pair<unsigned, unsigned> endPos,
							unsigned q);

bool insert_seeds_in_seedChain(seqan::String<seqan::Seed<seqan::Simple> > &seedChain,
										std::vector<std::pair<unsigned, unsigned> > &startPos,
										std::vector<std::pair<unsigned, unsigned> > &endPos,
										seqan::DnaString &seqOne, seqan::DnaString &seqTwo,
										unsigned q);

bool create_positionsContainers(std::vector<std::pair<unsigned, unsigned> > &startPos,
									std::vector<std::pair<unsigned, unsigned> > &endPos,
									seqan::String<seqan::Seed<seqan::Simple> > &seedChain,
									unsigned q, unsigned lenSeqOne, unsigned lenSeqTwo);

bool read_FASTA(seqan::CharString path, seqan::CharString  &id, seqan::DnaString &seq);

std::string get_file_contents(seqan::CharString filepath);

template<typename TInfix, typename TQueryId>
bool get_stellarMatches(seqan::StringSet<QueryMatches<StellarMatch<TInfix, TQueryId> > > & matches,
							seqan::SeedSet<seqan::Simple> &seedSet)
{
	typedef StellarMatch<TInfix, TQueryId> TMatch;
	typedef typename Size<typename TMatch::TAlign>::Type TSize;
	typedef typename Iterator<String<TMatch> >::Type TIterator;
	typedef Seed<Simple> SSeed;

	QueryMatches<TMatch> &queryMatches = value(matches, 0);

	TIterator it = begin(queryMatches.matches);
	TIterator itEnd = end(queryMatches.matches);

	int z = 0;
	if (it == itEnd)
	{
		std::cerr << "STELLAR: no matches found with entered options...." << std::endl;
		return false;
	}
	while (it != itEnd)
	{
		//debugging:
		std::cout << "beginPositionDatabase[" << z << "]: " << seqan::beginPosition((*it).row1) << std::endl;
		std::cout << "endPositionDatabase[" << z << "]: " << seqan::endPosition((*it).row1) << std::endl;
		std::cout << "beginPositionQuery[" << z << "]: " << seqan::beginPosition((*it).row2) << std::endl;
		std::cout << "endPositionQuery[" << z << "]: " << seqan::endPosition((*it).row2) << std::endl;

		SSeed stellarSeed(seqan::beginPosition((*it).row1),
				seqan::beginPosition((*it).row2),
				seqan::endPosition((*it).row1),
				seqan::endPosition((*it).row2));

		addSeed(seedSet, stellarSeed, seqan::Single());
		++it;
		++z;
	}
	return true;

}

bool write_seed_positions(std::vector< std::pair<unsigned, unsigned> > &array1,
							std::vector< std::pair<unsigned, unsigned> > &array2,
							const char* filepath);

bool output_alignment(seqan::CharString path,
						seqan::CharString &idOne, seqan::CharString &idTwo,
						seqan::Align<seqan::DnaString, seqan::ArrayGaps> &alignment,
						int score);

#endif  // #ifndef SANDBOX_MY_SANDBOX_APPS_LAGAN_NEU_LAGAN_FUNCTIONS_H_
