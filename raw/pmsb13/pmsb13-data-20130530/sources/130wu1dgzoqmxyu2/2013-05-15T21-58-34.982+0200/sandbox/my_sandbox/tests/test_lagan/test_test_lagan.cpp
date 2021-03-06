// ==========================================================================
//                                 test_lagan
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
// Author: Your Name <moritz.k@fu-berlin.de>
// ==========================================================================

#include <seqan/basic.h>
#include <seqan/file.h>
#include "test_test_lagan.h"
#include  "../../apps/lagan_omp/lagan_functions.h"

using namespace seqan;

SEQAN_DEFINE_TEST(test_my_app_funcs_create_positionsContainers)
{	
	unsigned q = 12;
	unsigned lenSeqOne = 200000;
	unsigned lenSeqTwo = 205000;
	
	typedef seqan::Seed<seqan::Simple> SSeed;
	seqan::String<SSeed> seedChain;
	SSeed seed(0u, 15u, 200u, 215u);
	SSeed seed2(50000u, 55000u, lenSeqOne - 20 * q, lenSeqTwo - 20 * q);	
	
	seqan::append(seedChain, seed);
	seqan::append(seedChain, seed2); 
	
	std::vector<std::pair<unsigned, unsigned> > startPositions;
	std::vector<std::pair<unsigned, unsigned> > endPositions;
	
	SEQAN_ASSERT_EQ(create_positionsContainers(startPositions,
								endPositions,
								seedChain,
								q, lenSeqOne, lenSeqTwo), 	true);
	
	SEQAN_ASSERT_EQ(startPositions.front().first, 200u);
	SEQAN_ASSERT_EQ(startPositions.front().second, 215u);
	SEQAN_ASSERT_EQ(endPositions.front().first, 50000u);
	SEQAN_ASSERT_EQ(endPositions.front().second, 55000u);
	
	SEQAN_ASSERT_EQ(startPositions.back().first, lenSeqOne - 20 * q);
	SEQAN_ASSERT_EQ(startPositions.back().second, lenSeqTwo - 20 * q);
	SEQAN_ASSERT_EQ(endPositions.back().first, lenSeqOne);
	SEQAN_ASSERT_EQ(endPositions.back().second, lenSeqTwo);
	
	clear(seedChain);
	
	SEQAN_ASSERT_EQ(create_positionsContainers(startPositions,
								endPositions,
								seedChain,
								q, lenSeqOne, lenSeqTwo), 	false);
								
	SSeed seed3(lenSeqOne - 15u, lenSeqTwo - 15u, lenSeqOne, lenSeqTwo);
	
	append(seedChain, seed3);
	startPositions.resize(0);
	endPositions.resize(0);
	
	SEQAN_ASSERT_EQ(create_positionsContainers(startPositions,
								endPositions,
								seedChain,
								q, lenSeqOne, lenSeqTwo), 	true);
								
	SEQAN_ASSERT_EQ(startPositions.front().first, 0u);
	SEQAN_ASSERT_EQ(startPositions.front().second, 0u);
	SEQAN_ASSERT_EQ(endPositions.front().first, lenSeqOne - 15u);
	SEQAN_ASSERT_EQ(endPositions.front().second, lenSeqTwo - 15u);
}

SEQAN_DEFINE_TEST(test_my_app_funcs_sortFunction)
{
	typedef seqan::Seed<seqan::Simple> SSeed;
	
	SSeed seed1(10, 50, 60, 100);
	SSeed seed2(20, 80, 100, 160);
	
	SEQAN_ASSERT_EQ(sortFunction(seed1, seed2), true);
	SEQAN_ASSERT_EQ(sortFunction(seed2, seed1), false);
}_

SEQAN_DEFINE_TEST(test_my_app_funcs_get_qGramArray)
{
	unsigned maxQ = 14;
	unsigned minQ = 4;
	unsigned parts = 1;
	
	std::vector<unsigned> qGramArray;
	
	qGramArray = get_qGramArray(maxQ, minQ, parts);
	
	SEQAN_ASSERT_EQ(qGramArray.size(), 2)
}
SEQAN_BEGIN_TESTSUITE(test_test_lagan)
{
    // Call tests.
	SEQAN_CALL_TEST(test_my_app_funcs_create_positionsContainers);
	SEQAN_CALL_TEST(test_my_app_funcs_sortFunction);
}
SEQAN_END_TESTSUITE
