// ==========================================================================
//                             adapterTrimming.h
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

#ifndef SANDBOX_GROUP3_APPS_SEQDPT_ADAPTERTRIMMING_H_
#define SANDBOX_GROUP3_APPS_SEQDPT_ADAPTERTRIMMING_H_

#include <seqan/align.h>

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// Align and strip paired-end reads.

template <typename TSeq>
seqan::Align<TSeq> alignPair(TSeq& seq1, TSeq& seq2)
{
	// When aligning the two sequences, the complementary sequence is reversed and
	// complemented, so we have an overlap alignment with complementary bases being the same.
	seqan::Dna5StringReverseComplement mod(seq2);
	typedef seqan::Align<TSeq> TAlign;
	TAlign align;
	seqan::resize(rows(align), 2);
	seqan::assignSource(row(align, 0), seq1);
	seqan::assignSource(row(align, 1), mod);

	// Overlap alignment. We again don't allow gaps, since they are unlikely in Illumina data.
	// The align configuration allows begin/end gaps. (On both sequences, different from adapter overlap)
    int score = globalAlignment(align, seqan::Score<int,seqan::Simple>(1,-1,-100),
    							seqan::AlignConfig<true, true, true, true>());
	std::cout << align << "\n";

	return align;
}

template <typename TRow>
int countTotalGaps(TRow& row)
{
	return length(row) - length(source(row));
}

template <typename TRow>
int getInsertSize(TRow& row1, TRow& row2)
{
	int seq1_length = length(source(row1));
	int seq2_length = length(source(row2));

	// Calculate overlap and overhangs.
	int overlap = seq1_length - countTotalGaps(row2);
	int seq2l = seqan::countGaps(seqan::begin(row1)); // Overhang of sequence 2 = Gaps at start of row 1.
	int seq1r = countTotalGaps(row2) - seqan::countGaps(seqan::begin(row2)); // Overhang of sequence 1 = Gaps at end of row 2.

	// Insert size: Add sequence lengths, subtract common region (overlap)
	// and subtract overhangs left and right of insert.
	return seq1_length + seq2_length - overlap - (seq2l + seq1r);
}

template <typename TSeq>
int stripPair(TSeq& seq1, TSeq& seq2)
{
	typedef seqan::Align<TSeq> TAlign;
	typedef typename seqan::Row<TAlign>::Type TRow;
	TAlign align = alignPair(seq1, seq2);
	TRow &row1 = row(align,0);
	TRow &row2 = row(align,1);

	// Use the overlap of the two sequences to determine the end position.
	int overlap = length(seq1) - countTotalGaps(row2);
	int mismatches = 0;

	// Get actual size of the insert (possible to determine from overlap).
	int insert_size = getInsertSize(row1, row2);
	std::cout << overlap << " " << insert_size << "\n";

	// If the sequences do not overlap (= 0) or the overlap is not big enough
	// to make one of the reads go into an adapter sequence -> do not cut.
	if (overlap < insert_size)
		return insert_size;

	// Now cut both sequences to insert size (no cuts happen if they are smaller)
	seqan::erase(seq1, insert_size, length(seq1));
	seqan::erase(seq2, insert_size, length(seq2));

	return insert_size;
}

// Align and strip single reads.

template <typename TSeq, typename TAdapter>
seqan::Pair<int, int> alignAdapter(TSeq& seq, TAdapter& adapter)
{
	typedef seqan::Align<TSeq> TAlign;
	typedef typename seqan::Row<TAlign>::Type TRow;
	TAlign align;
	seqan::resize(rows(align), 2);
	seqan::assignSource(row(align, 0), seq);
	seqan::assignSource(row(align, 1), adapter);

	TRow &row1 = row(align,0);
	TRow &row2 = row(align,1);

    // Global free-end alignment. The Alignment configuration specifies that gaps on the
    // start of the read and at the end of the adapter are not penalized -> 5'-3' overlap.
    // We also don't allow gaps by setting the gap penalty high.
    int score = globalAlignment(align, seqan::Score<int,seqan::Simple>(1,-1,-100),
    							seqan::AlignConfig<true, false, true, false>(), 0, length(seq));

    //Overlap is sequence length - gaps before the adapter, since we don't allow any gaps later on.
    int overlap = length(seq)-seqan::countGaps(begin(row2));
    int mismatches = (overlap-score)/2;

	return seqan::Pair<int,int>(overlap, mismatches);
}

template <typename TSeq>
seqan::Pair<int, int> alignRevComplAdapter(TSeq& seq, TSeq& adapter)
{
	seqan::Dna5StringReverseComplement mod(adapter);
	return alignAdapter(seq, mod);
}

// Manually set error rate and minimum length.
template <typename TSeq>
int stripAdapter(TSeq& seq, TSeq& adapter, int min_length, int errors)
{
	seqan::Pair<int,int> res = alignAdapter(seq, adapter);
	if (res.i1 >= min_length && res.i2 <= errors)
	{
		seqan::erase(seq, length(seq) - res.i1, length(seq));
		return res.i1;
	}

	return 0;
}

// Auto configuration.
template <typename TSeq>
int stripAdapter(TSeq& seq, TSeq& adapter)
{
	seqan::Pair<int,int> res = alignAdapter(seq, adapter);
	int errors = 0; // No errors for overlap up to 5 bases.
	if (res.i1 > 5)
		errors = 1; // One error for overlap up to 10 bases.
	if (res.i1 > 10)
		errors = 0.15*res.i1; //15% of overlap length otherwise.

	// If there are fewer mismatches than allowed -> cut.
	if (res.i2 <= errors)
	{
		seqan::erase(seq, length(seq) - res.i1, length(seq));
		return res.i1;
	}

	return 0;
}

#endif  // #ifndef SANDBOX_GROUP3_APPS_SEQDPT_ADAPTERTRIMMING_H_
