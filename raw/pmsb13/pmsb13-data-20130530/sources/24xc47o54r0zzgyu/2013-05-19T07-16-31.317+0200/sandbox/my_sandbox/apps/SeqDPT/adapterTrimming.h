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
// Author: Benjamin Strauch
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

// Custom scoring matrix for ignoring alignments with N.
namespace seqan{
	/*!
	 * \brief Struct used to define a new custom scoring matrix.
	 */
	struct AdapterScoringMatrix {};

	/*!
	 * \brief Struct containing data for the ::AdapterScoringMatrix custom scoring matrix.
	 * Matches score 1, mismatches -1 and matches against N with 0.
	 */
	template <>
	struct ScoringMatrixData_<int, Dna5, AdapterScoringMatrix> {
		enum {
			VALUE_SIZE = ValueSize<Dna5>::VALUE,
			TAB_SIZE = VALUE_SIZE * VALUE_SIZE
		};

		static inline int const * getData() {
			static int const _data[TAB_SIZE] = {
			   1, -1, -1, -1, 0,
			  -1,  1, -1, -1, 0,
			  -1, -1,  1, -1, 0,
			  -1, -1, -1,  1, 0,
			   0,  0,  0,  0, 0
			};
			return _data;
		}
	};
}

/*!
 * \brief A struct encapsulating information about the match algorithm.
 */
struct Mode {};
/*!
 * \brief Tagging struct for the automatic match algorithm.
 */
struct Auto : Mode {};
/*!
 * \brief Tagging struct representing the the match algorithm working
 * with values supplied by the user. Saves those values as members.
 */
struct User : Mode {
	int min_length; /*!< The minimum length of the overlap. */
	int errors;     /*!< The maximum number of errors we allow. */
	User (int m, int e): min_length(m), errors(e) {}
};

// ============================================================================
// Metafunctions
// ============================================================================

/*!
 * \brief A metafunction which constructs a ModifiedString type for an Alphabet.
 */
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

/*!
 * \brief Align two sequences, returning the score and align object of the alignment.
 * \param seq1 The first sequence of the alignment.
 * \param seq2 The second sequence of the alignment.
 * \param config The alignment configuration used by the Seqan's globalAlignment method.
 * \param band Whether the alignment's lower diagonal should be banded at -2. (Two diagonals
 * below the main diagonal.) Useful to speed up computation for 5'-3' overlap alignments.
 * \remark The main intended uses are AlignConfig<true, false, true, false> with band = true
 * and AlignConfig<true, true, true, true> with band = false.
 * \remark We use a custom scoring scheme which scores 1 for match, -1 for mismatch, 0 for match with N.
 * \return A pair containing the score of the alignment and the aling object.
 */
template <typename TSeq1, typename TSeq2, bool TTop, bool TLeft, bool TRight, bool TBottom>
seqan::Pair<unsigned, seqan::Align<TSeq1> > alignPair(TSeq1& seq1, TSeq2& seq2,
		const seqan::AlignConfig<TTop, TLeft, TRight, TBottom>& config, bool band = false)
{
	typedef seqan::Align<TSeq1> TAlign;
	TAlign align;
	seqan::resize(rows(align), 2);
	seqan::assignSource(row(align, 0), seq1);
	seqan::assignSource(row(align, 1), seq2);

	// Overlap alignment. We again don't allow gaps, since they are unlikely in Illumina data.
	// The align configuration allows begin/end gaps. (On both sequences, different from adapter overlap)
	Score<int, ScoreMatrix<Dna5, AdapterScoringMatrix> > adapterScore;
	unsigned score;
    if (band)
    	score = globalAlignment(align, adapterScore, config, -2, length(seq1));
    else
    	score = globalAlignment(align, adapterScore, config);

	return seqan::Pair<unsigned, TAlign>(score, align);
}

/*!
 * \brief Returns the total number of gaps in a row object.
 * \param row The gapped row object used by seqan alignments.
 */
template <typename TRow>
unsigned countTotalGaps(TRow& row)
{
	return length(row) - length(source(row));
}

/*!
 * \brief Determines the overlap between two (free-shift) aligned gapless sequences.
 * \param align An alignment object with two aligned overlapping sequences.
 * \pre The aligned sequences must not contain any internal gaps.
 * \return The number of overlapping positions.
 */
template <typename TAlign>
unsigned getOverlap(TAlign& align)
{
	typedef typename seqan::Row<TAlign>::Type TRow;
	TRow &row1 = seqan::row(align,0);
	TRow &row2 = seqan::row(align,1);

	return length(source(row1)) - countTotalGaps(row2);
}

/*!
 * \brief Given an alignment with two overlapping (forward and reverse) sequences,
 * this function determines the actual size of the insert they cover.
 * \param align An alignment object with two aligned overlapping sequences.
 * \pre The aligned sequences must overlap and not contain any internal gaps.
 * \remark This method can reconstruct the insert perfectly, provided the alignment
 * object represents the real alignment. There can be small discrepancies if indels occurred.
 * \return The insert size that was determined from the overlap alignment.
 */
template <typename TAlign>
unsigned getInsertSize(TAlign& align)
{
	typedef typename seqan::Row<TAlign>::Type TRow;
	TRow &row1 = seqan::row(align,0);
	TRow &row2 = seqan::row(align,1);

	unsigned seq1_length = length(source(row1));
	unsigned seq2_length = length(source(row2));

	// Calculate overlap and overhangs.
	unsigned overlap = seq1_length - countTotalGaps(row2);
	unsigned seq2l = seqan::countGaps(seqan::begin(row1)); // Overhang of sequence 2 = Gaps at start of row 1.
	unsigned seq1r = countTotalGaps(row2) - seqan::countGaps(seqan::begin(row2)); // Overhang of sequence 1 = Gaps at end of row 2.

	// Insert size: Add sequence lengths, subtract common region (overlap)
	// and subtract overhangs left and right of insert.
	return seq1_length + seq2_length - overlap - (seq2l + seq1r);
}

/*!
 * \brief Removes adapter contamination from paired-end reads.
 * \param seq1 The forward read.
 * \param seq2 The backward read.
 * \remark This method is very accurate and does not need any knowledge of
 * the specific adapter sequences contaminating the reads.
 * \return The determined actual insert size or 0 if the sequences don't overlap.
 * \remark The number of trimmed bases is len(seq) - insert, if the sequence was greater than the insert size.
 */
template <typename TSeq>
unsigned stripPair(TSeq& seq1, TSeq& seq2)
{
	// When aligning the two sequences, the complementary sequence is reversed and
	// complemented, so we have an overlap alignment with complementary bases being the same.
	typedef typename Value<TSeq>::Type TAlphabet;
	typedef typename STRING_REVERSE_COMPLEMENT<TAlphabet>::Type TReverseComplement;
	TReverseComplement mod(seq2);

	typedef seqan::Align<TSeq> TAlign;
	seqan::Pair<unsigned, TAlign> ret = alignPair(seq1, mod, seqan::AlignConfig<true, true, true, true>());
	unsigned score = ret.i1;
	TAlign align = ret.i2;

	// Use the overlap of the two sequences to determine the end position.
	unsigned overlap = getOverlap(align);
	unsigned mismatches = (overlap-score)/2;

	// We require a certain correct overlap to exclude spurious hits.
	// (Especially reverse 3'->5' alignments not caused by adapters.)
	if (overlap <= 5 || mismatches > overlap*0.15)
		return 0;

	// Get actual size of the insert (possible to determine from overlap etc.).
	unsigned insert = getInsertSize(align);

	// Now cut both sequences to insert size (no cuts happen if they are smaller)
	if (length(seq1) > insert)
		seqan::erase(seq1, insert, length(seq1));
	if (length(seq2) > insert)
		seqan::erase(seq2, insert, length(seq2));

	return insert;
}

/*!
 * \brief Removes adapter contamination from paired-end reads with adapter information.
 * \param seq1 The forward read.
 * \param adapter1 The adapter that contaminates the forward read.
 * \param seq2 The backward read.
 * \param adapter2 The adapter whose reverse complement contaminates the backward read.
 * \remark This method can be less accurate than the overload which doesn't use any adapter
 * sequence information, since it is more constrained in how it tries to overlap the sequences.
 * \return The determined actual insert size or 0 if the sequences don't overlap.
 * \remark The number of trimmed bases is len(seq) - insert, if the sequence was greater than the insert size.
 */
template <typename TSeq>
unsigned stripPair(TSeq& seq1, TSeq& adapter1, TSeq& seq2, TSeq& adapter2)
{
	typedef typename Value<TSeq>::Type TAlphabet;
	typedef typename STRING_REVERSE_COMPLEMENT<TAlphabet>::Type TReverseComplement;

	// Add adapters to the sequences.
	seqan::insert(seq1, 0, adapter2);
	seqan::insert(seq2, 0, adapter1);

	TReverseComplement mod(seq2);

	typedef seqan::Align<TSeq> TAlign;
	seqan::Pair<unsigned, TAlign> ret = alignPair(seq1, mod, seqan::AlignConfig<true, false, true, false>(), true);
	TAlign align = ret.i2;
	unsigned score = ret.i1;
	unsigned overlap = getOverlap(align);
	unsigned mismatches = (overlap-score)/2;

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
	if (mismatches >= overlap*0.15)
		return 0;

	// Cut down to insert size, if necessary.
	if (length(seq1) > insert)
		seqan::erase(seq1, insert, length(seq1));
	if (length(seq2) > insert)
		seqan::erase(seq2, insert, length(seq2));

	return insert;
}

/*!
 * \brief Removes adapter contamination from paired-end reads.
 * \param seq1 The set of forward reads.
 * \param seq2 The set of backward reads.
 * \remark This method is very accurate and does not need any knowledge of
 * the specific adapter sequences contaminating the reads.
 * \remark This method can operate concurrently with OpenMP.
 * \return The number of sequences that had some bases removed.
 */
template <typename TSeq>
unsigned stripPairBatch(seqan::StringSet<TSeq>& set1, seqan::StringSet<TSeq>& set2)
{
	unsigned stripped = 0;
	int len = length(set1);
	#pragma omp parallel for reduction(+:stripped)
	for (int i=0; i < len; ++i)
	{
		if (stripPair(value(set1, i), value(set2, i)))
			++stripped;
	}

	return stripped;
}

/*!
 * \brief Aligns a sequence to an adapter.
 * \param seq The sequence whose adapter contamination should be removed.
 * \param adapter The adapter that might contaminate the sequence.
 * \return A pair containing the score of the alignment and the alignment object.
 */
template <typename TSeq, typename TAdapter>
seqan::Pair<unsigned, seqan::Align<TSeq> > alignAdapter(TSeq& seq, TAdapter& adapter)
{
	// Global free-end alignment. The Alignment configuration specifies that gaps on the
	// start of the read and at the end of the adapter are not penalized -> 5'-3' overlap.
	// We also don't allow gaps by setting the gap penalty high.
	return alignPair(seq, adapter, seqan::AlignConfig<true, false, true, false>(), true);
}

/*!
 * \brief Aligns a sequence to the reverse complement of an adapter.
 * \param seq The sequence whose adapter contamination should be removed.
 * \param adapter The adapter whose reverse complement might contaminate the sequence.
 * \return A pair containing the score of the alignment and the alignment object.
 */
template <typename TSeq>
seqan::Pair<unsigned, seqan::Align<TSeq> > alignRevComplAdapter(TSeq& seq, TSeq& adapter)
{
	typedef typename Value<TSeq>::Type TAlphabet;
	typedef typename STRING_REVERSE_COMPLEMENT<TAlphabet>::Type TReverseComplement;
	TReverseComplement mod(adapter);

	return alignAdapter(seq, mod);
}

/*!
 * \brief Checks if a overlap of an alignment is accepted,
 * based on mismatches and the length of the overlap.
 * \param overlap The number of overlapping positions in the overlap alignment.
 * \param mismatches The number of allowed mismatches in the overlapping region.
 * \param Test
 * \remark This method automatically uses a very simple heuristic to determine matches.
 * \return Bool indicating if the alignment is significant.
 */
bool isMatch(int overlap, int mismatches, const Auto &)
{
	int errors = 0; // No errors for overlap up to 5 bases.
	if (overlap > 5)
		errors = 1; // One error for overlap up to 10 bases.
	if (overlap > 10)
		errors = 0.33*overlap; //15% of overlap length otherwise.

	return mismatches <= errors;
}

/*!
 * \brief Checks if a overlap of an alignment is accepted,
 * based on mismatches and the length of the overlap.
 * \param overlap The number of overlapping positions in the overlap alignment.
 * \param mismatches The number of allowed mismatches in the overlapping region.
 * \param userOptions Parameters specifying requirements for the overlap.
 * \remark This overload is used to let the user specify match requirements.
 * \return Bool indicating if the alignment is significant.
 */
bool isMatch(int overlap, int mismatches, const User& userOptions)
{
	return overlap >= userOptions.min_length && mismatches <= userOptions.errors;
}

/*!
 * \brief Remove adapter sequence from a sequence.
 * \param seq The sequence whose adapter contamination should be removed.
 * \param adapter The adapter sequence that might contaminate the sequence.
 * \param spec The match algorithm that decides whether a match was significant.
 * \return The overlap of the sequence with the adapter.
 */
template <typename TSeq, typename TAdapter, typename TSpec>
unsigned stripAdapter(TSeq& seq, TAdapter& adapter, TSpec& spec)
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


/*!
 * \brief Remove adapter sequence from a set of sequences.
 * \param set A StringSet of sequences whose adapter contaminations should be removed.
 * \param adapter The adapter sequence that might contaminate the sequences.
 * \param spec The match algorithm that decides whether a match was significant.
 * \remark This method simply applies ::stripAdapter to all sequences in ::set.
 * \remark This method can operate concurrently with OpenMP.
 * \return The number of sequences that had some bases removed.
 */
template <typename TSeq, typename TAdapter, typename TSpec>
unsigned stripAdapterBatch(seqan::StringSet<TSeq>& set, TAdapter& adapter, TSpec const & spec)
{
	unsigned stripped = 0;
	int len = length(set);
	#pragma omp parallel for reduction(+:stripped)
	for (int i=0; i < len; ++i)
	{
		if (stripAdapter(value(set, i), adapter, spec))
			++stripped;
	}

	return stripped;
}

/*!
 * \brief Simple interface to align the reverse complement of an adapter to a batch of sequences.
 * to a set of reads and remove significant matches.
 * \param set A StringSet of sequences whose adapter contaminations should be removed.
 * \param adapter The adapter sequence whose reverse complement might contaminate the sequences.
 * \param spec The match algorithm that decides whether a match was significant.
 * \remark This interface can be used to detect contamination in reverse reads of paired-end
 * 	reads. Those reads might read into the reverse complement of an adapter.
 * \return The number of sequences that had some bases removed.
 */
template <typename TSeq, typename TSpec>
unsigned stripReverseAdapterBatch(seqan::StringSet<TSeq>& set, TSeq& adapter, TSpec const & spec)
{
	typedef typename Value<TSeq>::Type TAlphabet;
	typedef typename STRING_REVERSE_COMPLEMENT<TAlphabet>::Type TReverseComplement;
	TReverseComplement mod(adapter);

	return stripAdapterBatch(set, mod, spec);
}

#endif  // #ifndef SANDBOX_GROUP3_APPS_SEQDPT_ADAPTERTRIMMING_H_
