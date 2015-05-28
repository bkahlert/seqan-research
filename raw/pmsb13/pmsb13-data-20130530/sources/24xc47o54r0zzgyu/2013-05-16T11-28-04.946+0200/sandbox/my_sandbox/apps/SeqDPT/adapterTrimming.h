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
#include <seqan/find.h>

#define SCORING seqan::SimpleScore(1,-1,-100)
using namespace seqan;

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// Adapter trimming configuration.
struct Mode {};
struct Auto : Mode {};
struct User : Mode {
	int min_length;
	int errors;
	User (int m, int e): min_length(m), errors(e) {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// Returns the appropriate Modified String for a reverse complement view.
template <class TValue>
struct STRING_REVERSE_COMPLEMENT
{
	typedef ModifiedString<
			ModifiedString<	String<TValue>, ModView< FunctorComplement<TValue> > >,
			ModReverse>	Type;
};

// ============================================================================
// Functions
// ============================================================================

/**
 * Simply align two sequences with the specified alignment configuration.
 * A constant scoring scheme is used.
 */
template <typename TSeq1, typename TSeq2, bool TTop, bool TLeft, bool TRight, bool TBottom>
seqan::Pair<int, seqan::Align<TSeq1> > alignPair(TSeq1& seq1, TSeq2& seq2, seqan::AlignConfig<TTop, TLeft, TRight, TBottom> config)
{
	typedef seqan::Align<TSeq1> TAlign;
	TAlign align;
	seqan::resize(rows(align), 2);
	seqan::assignSource(row(align, 0), seq1);
	seqan::assignSource(row(align, 1), seq2);

	// Overlap alignment. We again don't allow gaps, since they are unlikely in Illumina data.
	// The align configuration allows begin/end gaps. (On both sequences, different from adapter overlap)
    int score = globalAlignment(align, SCORING, config);
	//std::cout << align << "\n";

	return seqan::Pair<int, TAlign>(score, align);
}

/**
 * Gets the total number of gaps in the alignment row as the difference
 * of the original sequence and the larger (because of the gaps) aligned sequence.
 */
template <typename TRow>
int countTotalGaps(TRow& row)
{
	return length(row) - length(source(row));
}

/**
 * Determines overlap between the two alignment rows, assuming
 * that the alignment only contains leading and trailing gaps.
 */
template <typename TAlign>
int getOverlap(TAlign& align)
{
	typedef typename seqan::Row<TAlign>::Type TRow;
	TRow &row1 = seqan::row(align,0);
	TRow &row2 = seqan::row(align,1);

	return length(source(row1)) - countTotalGaps(row2);
}

/**
 * Determines actual insert size. Only works if sequences actually overlap.
 */
template <typename TAlign>
int getInsertSize(TAlign& align)
{
	typedef typename seqan::Row<TAlign>::Type TRow;
	TRow &row1 = seqan::row(align,0);
	TRow &row2 = seqan::row(align,1);

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

/**
 * Strips adapter sequences from paired reads. Sequence 1 is aligned to the
 * reverse complementary version of Sequence 2. If they overlap, we can recover
 * the original insert size and cut both sequences to that size, since anything
 * beyond that insert size must be some kind of adapter sequence.
 */
template <typename TSeq>
int stripPair(TSeq& seq1, TSeq& seq2)
{
	// When aligning the two sequences, the complementary sequence is reversed and
	// complemented, so we have an overlap alignment with complementary bases being the same.
	STRING_REVERSE_COMPLEMENT<seqan::Dna5Q >::Type mod(seq2);

	typedef seqan::Align<TSeq> TAlign;
	seqan::Pair<int, TAlign> ret = alignPair(seq1, mod, seqan::AlignConfig<true, true, true, true>());
	int score = ret.i1;
	TAlign align = ret.i2;

	// Use the overlap of the two sequences to determine the end position.
	int overlap = getOverlap(align);
	int mismatches = (overlap-score)/2;

	// We require a certain correct overlap to exclude spurious hits.
	// (Especially reverse 3'->5' alignments not caused by adapters.)
	if (overlap <= 5 || mismatches > overlap*0.3)
		return 0;

	// Get actual size of the insert (possible to determine from overlap etc.).
	int insert = getInsertSize(align);

	// Now cut both sequences to insert size (no cuts happen if they are smaller)
	if (length(seq1) > insert)
		seqan::erase(seq1, insert, length(seq1));
	if (length(seq2) > insert)
		seqan::erase(seq2, insert, length(seq2));

	return insert;
}

/**
 * Pairwise trimming with adapter information.
 */
template <typename TSeq>
int stripPair(TSeq& seq1, TSeq& adapter1, TSeq& seq2, TSeq& adapter2)
{
	typedef STRING_REVERSE_COMPLEMENT<seqan::Dna5Q>::Type modType;

	//std::cout << seq1 << "\n" << adapter2 << "\n-----\n" << seq2 << "\n" << adapter1 << "\n~~~~~~~~~~\n";

	// Add adapters to the sequences.
	seqan::insert(seq1, 0, adapter2);
	seqan::insert(seq2, 0, adapter1);

	modType mod(seq2);

	typedef seqan::Align<TSeq> TAlign;
	seqan::Pair<int, TAlign> ret = alignPair(seq1, mod, seqan::AlignConfig<true, false, true, false>());
	TAlign align = ret.i2;
	int score = ret.i1;
	int overlap = getOverlap(align);
	int mismatches = (overlap-score)/2;

	// We can still get the proper insert size from the overlap, since we only
	// get properly overlapping alignments in this alignment mode. We just have
	// to subtract the adapter additions. One exception is the case where they
	// don't overlap at all, but in this case getInsertSize returns 0.
	int insert = getInsertSize(align);
	if (insert > 0)
		insert = insert - length(adapter1) - length(adapter2);

	// Remove adapter sequences again.
	seqan::erase(seq1, 0, length(adapter2));
	seqan::erase(seq2, 0, length(adapter1));

	// Check the quality of the overlap.
	if (mismatches >= overlap*0.3)
		return 0;

	// Cut down to insert size, if necessary.
	if (length(seq1) > insert)
		seqan::erase(seq1, insert, length(seq1));
	if (length(seq2) > insert)
		seqan::erase(seq2, insert, length(seq2));

	return insert;
}

/**
 * Remove adapters from a set of paired-end reads.
 */
template <typename TSeq>
int stripPairBatch(seqan::StringSet<TSeq>& set1, seqan::StringSet<TSeq>& set2)
{
	int stripped = 0;
	for (int i=0; i < length(set1); ++i)
	{
		if (stripPair(value(set1, i), value(set2, i))) ++stripped;
	}

	return stripped;
}

/**
 * Align a sequence to its adapter. This is done by doing an overlap alignment.
 */
template <typename TSeq, typename TAdapter>
seqan::Pair<int, seqan::Align<TSeq> > alignAdapter(TSeq& seq, TAdapter& adapter)
{
	// Global free-end alignment. The Alignment configuration specifies that gaps on the
	// start of the read and at the end of the adapter are not penalized -> 5'-3' overlap.
	// We also don't allow gaps by setting the gap penalty high.
	return alignPair(seq, adapter, seqan::AlignConfig<true, false, true, false>());
}

/**
 * Align a sequence to a reversed adapter. This is done by doing an overlap
 * alignment with a reversed and complemented view of the adapter sequence.
 */
template <typename TSeq>
seqan::Pair<int, seqan::Align<TSeq> > alignRevComplAdapter(TSeq& seq, TSeq& adapter)
{
	STRING_REVERSE_COMPLEMENT<Value<TSeq> > mod(adapter);
	return alignAdapter(seq, mod);
}

/**
 * Check if the overlap alignment is accepted, based on
 * mismatches and the length of the overlap.
 *
 * Automatic configuration.
 */
bool isMatch(int overlap, int mismatches, Auto)
{
	int errors = 0; // No errors for overlap up to 5 bases.
	if (overlap > 5)
		errors = 1; // One error for overlap up to 10 bases.
	if (overlap > 10)
		errors = 0.15*overlap; //15% of overlap length otherwise.

	return mismatches <= errors;
}

/**
 * Check if the overlap alignment is accepted, based on
 * mismatches and the length of the overlap.
 *
 * Values for minimum length and maximum errors given by the user.
 */
bool isMatch(int overlap, int mismatches, User u)
{
	return overlap >= u.min_length && mismatches <= u.errors;
}

/**
 * Remove adapter sequence from a sequence.
 */
template <typename TSeq, typename TSpec>
int stripAdapter(TSeq& seq, TSeq& adapter, TSpec& spec)
{
	typedef seqan::Align<TSeq> TAlign;
	seqan::Pair<int, TAlign> ret = alignAdapter(seq, adapter);
	int overlap = getOverlap(ret.i2);
	int score = ret.i1;
	int mismatches = (overlap-score)/2;

	if (isMatch(overlap, mismatches, spec))
	{
		seqan::erase(seq, length(seq) - overlap, length(seq));
		return overlap;
	}

	return 0;
}

/**
 * Remove adapter sequence from a set of sequences.
 */
template <typename TSeq, typename TSpec>
int stripAdapterBatch(seqan::StringSet<TSeq>& set, TSeq& adapter, TSpec const & spec)
{
	int stripped = 0;
	for (int i=0; i < length(set); ++i)
	{
		if (stripAdapter(value(set, i), adapter, spec)) ++stripped;
	}

	return stripped;
}

#endif  // #ifndef SANDBOX_GROUP3_APPS_SEQDPT_ADAPTERTRIMMING_H_
