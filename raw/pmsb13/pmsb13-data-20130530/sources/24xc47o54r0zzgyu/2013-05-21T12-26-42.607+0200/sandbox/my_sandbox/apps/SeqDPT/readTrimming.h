// ==========================================================================
//                               readTrimming.h
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

/*! \file readTrimming.h
\brief Contains the functions for read trimming.
*/

#ifndef SANDBOX_GROUP3_APPS_SEQDPT_READTRIMMING_H_
#define SANDBOX_GROUP3_APPS_SEQDPT_READTRIMMING_H_

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// Tags for choosing the appropriate trimming method.
/*!
 * \brief The tagging structure for the Tail trimming algorithm.
 */
struct Tail {};

/*!
 * \brief The tagging structure for the BWA trimming algorithm.
 */
struct BWA {};

/*!
 * \brief The tagging structure for the window trimming algorithm.
 */
struct Mean {
	unsigned window;	/*!< The window size used for quality calculations. */
	Mean(unsigned w) : window(w) {}
};

/*!
 * \brief This structure wraps a sequence with dedicated Dna5 and quality
 * strings so we can use it with our trimming function.
 */
struct Dna5QAdapter{
	seqan::String<seqan::Dna5>& seq; /*!< The Dna5 string of the sequence. */
	seqan::CharString& qual; 		 /*!< The quality string of the sequence. */
	Dna5QAdapter(seqan::String<seqan::Dna5> &s, seqan::CharString &q): seq(s), qual(q) {}
};

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

/*!
 * \brief Determines the quality at position i of the Dna5QString.
 * \param seq The Dna5QString containing bases and qualities.
 * \param i The index of the base whose quality shall be returned.
 * \return Phred quality of the base at position \b i.
 */
inline unsigned getQuality(const seqan::String<seqan::Dna5Q>& seq, unsigned i)
{
	return seqan::getQualityValue(seq[i]);
}

/*!
 * \brief Determines the quality at position i in the quality string
 * of the ::Dna5QAdapter.
 * \param seq The ::Dna5QAdapter structure containing the quality string.
 * \param i The index of the base whose quality shall be returned.
 * \return Phred quality of the base at position \b i.
 */
inline unsigned getQuality(const Dna5QAdapter& seq, unsigned i)
{
	return seqan::ordValue(seq.qual[i]) - 33;
}


/*!
 * \brief A specialization of the length function for use with the wrapper
 * ::Dna5QAdapter. Returns the length of the sequence contained in the wrapper.
 */
unsigned length(const Dna5QAdapter & seq)
{
	return length(seq.seq);
}

/*!
 * \brief A specialization of the erase function for use with the wrapper
 * ::Dna5QAdapter. Simply trimmes both sequence and quality strings separately.
 */
void erase(const Dna5QAdapter & seq, unsigned begin, unsigned end)
{
	seqan::erase(seq.seq , begin, end);
	seqan::erase(seq.qual, begin, end);
}

// Trimming methods
// ----------------------------------------------------------------------------
/*!
 * \brief A very simple trimming mechanism that removes low quality bases from the end.
 * \param seq The sequence that shall be trimmed.
 * \param cutoff The minimum quality required.
 * \return The trimming position, i.e. the first base to be removed.
 * \remark Simply cut off as many low quality bases from the end as possible before finding a good one.
 */
template <typename TSeq>
unsigned _trimRead(const TSeq& seq, unsigned const cutoff, Tail const &)
{
	unsigned cut_pos = length(seq);
	for (int i=length(seq)-1; i >= 0; --i)
		if (getQuality(seq, i) < cutoff)
			cut_pos = i;

	return cut_pos;
}

/*!
 * \brief A trimming algorithm using the method implemented in BWA.
 * \param seq The sequence that shall be trimmed.
 * \param cutoff The minimum quality required.
 * \return The trimming position, i.e. the first base to be removed.
 * \remark Trimming mechanism used in BWA. Trim to argmax_x sum_{i=x+1}^l {cutoff - q_i}
 */
template <typename TSeq>
unsigned _trimRead(const TSeq& seq, unsigned const cutoff, BWA const &)
{
	int max_arg = length(seq), sum = 0, max = 0;
	for (int i=length(seq)-1; i >= 0; --i)
	{
		sum += cutoff - getQuality(seq, i);
		if (sum < 0)
			break;
		if (sum > max)
		{
			max = sum;
			max_arg = i;
		}
	}

	return max_arg;
}

/*!
 * \brief A trimming algorithm using a sliding window method to find the trimming position.
 * \param seq The sequence that shall be trimmed.
 * \param cutQual The minimum quality required.
 * \param spec The algorithm tag containing the window size used for trimming.
 * \return An unsigned int representing the trimming position, i.e. the position of the first base to be removed.
 */
// Trim by shifting a window over the sequence and cut where avg. qual. in window turns bad the first time.
template <typename TSeq>
unsigned _trimRead(const TSeq& seq, unsigned const cutQual, Mean const & spec)
{
	unsigned window = spec.window;
	unsigned avg = 0, i = 0;
	// Work with absolute cutoff in window to avoid divisions.
	unsigned cutoff = cutQual*window;

	// Calculate average quality of initial window.
	for (i=0; i < window; ++i)
		avg += getQuality(seq, i);

	// Shift window over read and keep mean quality, update in constant time.
	for (i=0; i < length(seq) && avg >= cutoff; ++i)
	{
		// Take care only not to go over the end of the sequence. Shorten window near the end.
		avg -= getQuality(seq, i);
		avg += i + window < length(seq) ? getQuality(seq, i + window) : 0;
	}

	// i now holds the start of the first window that turned bad.
	return i;
}

/*!
 * \brief Interface that trims a sequence.
 * \param seq The sequence that shall be trimmed.
 * \param cutoff The minimum quality required from a base.
 * \param spec The trimming algorithm used for trimming.
 * \return The number of bases trimmed from the sequence.
 * \remark This function takes dedicated sequence and quality strings,
 * puts them into a wrapper structure and calls ::trimRead(TSeq& seq, unsigned const cutoff, TSpec const & spec).
 */
template <typename TSeq, typename TQual, typename TSpec>
unsigned trimRead(TSeq& seq, TQual& qual, unsigned const cutoff, TSpec const & spec)
{
	Dna5QAdapter a = Dna5QAdapter(seq, qual);
	return trimRead(a, cutoff, spec);
}

/*!
 * \brief Interface that trims a sequence.
 * \param seq The sequence that shall be trimmed.
 * \param cutoff The minimum quality required from a base.
 * \param spec The trimming algorithm used for trimming.
 * \return The number of bases trimmed from the sequence.
 */
template <typename TSeq, typename TSpec>
unsigned trimRead(TSeq& seq, unsigned const cutoff, TSpec const & spec)
{
	unsigned ret, cut_pos;
	cut_pos = _trimRead(seq, cutoff, spec);

	ret = length(seq) - cut_pos;
	erase(seq, cut_pos, length(seq));

	return ret;
}

/*!
 * \brief Trims a set of reads and marks reads that are too short for removal.
 * \param removedID A collection that stores the IDs of sequences marked for removal.
 * \param seqSet A collection of sequences that shall be trimmed.
 * \param cutoff The minimum quality required from a base.
 * \param min_length The minimum length required after trimming.
 * Shorter sequences are marked for removal in \b removedID.
 * \param spec The trimming algorithm used for trimming.
 * \return The number of reads where bases were removed.
 */
template <typename TRem, typename TSet, typename TSpec>
unsigned _trimReads(TRem & removedID, TSet & seqSet, unsigned const cutoff,
				  unsigned const min_length, TSpec const & spec)
{
	typedef typename seqan::Value<TSet>::Type TSeq;

	int trimmedReads = 0;
	int len = length(seqSet);

	#pragma omp parallel for reduction(+:trimmedReads)
	for (int i=0; i < len; ++i)
	{
		TSeq& read = value(seqSet, i);
		unsigned trimmed = trimRead(read, cutoff, spec);
		if (trimmed > 0) ++trimmedReads;
		if (length(read) < min_length){
			#pragma omp critical
			seqan::append(removedID, seqan::positionToId(seqSet, i));
		}
	}

	std::sort(seqan::begin(removedID), seqan::end(removedID));
	return trimmedReads;
}

/*!
 * \brief Removes two sequences at the position \b id from two StringSets.
 * \param set1 The first set where position \b id is removed.
 * \param set2 The second set where position \b id is removed.
 * \param id The id that identifies the read to be removed.
 */
template <typename TSet1, typename TSet2, typename TId>
void removeValues(TSet1 & set1, TSet2 & set2, TId id)
{
	seqan::removeValueById( set1, id);
	seqan::removeValueById( set2, id);
}

/*!
 * \brief Trims bad quality bases from a set of sequences.
 * \param idSet StringSet of FastA-IDs for the first set of sequences.
 * \param seqSet StringSet containing the forward reads of the paired sequences.
 * \param cutoff The minimum quality required of a base.
 * \param min_length The minimum length required after trimming. Shorter sequences
 *  will be deleted.
 * \param spec The trimming algorithm used to trim the sequences.
 * \return The number of sequences which had bases removed from.
 */
template <typename TId, typename TSeq, typename TSpec>
unsigned trimBatch(seqan::StringSet<TId> & idSet, seqan::StringSet<TSeq> & seqSet,
		unsigned const cutoff, unsigned const min_length, TSpec const & spec)
{
	// Only save the StringSet-Ids of the strings that are going to be
	// removed. This avoids index problems when removing from the StringSet
	// in a running loop and enables easy parallelization with OpenMP.
	typedef typename seqan::Id<seqan::StringSet<TSeq> >::Type TRem;
	seqan::String<TRem> removedID;
	unsigned trimmedReads = _trimReads(removedID, seqSet, cutoff, min_length, spec);

	// Remove deleted sequences from StringSet.
	for (int i=seqan::length(removedID)-1; i >= 0 ; --i)
		removeValues(idSet, seqSet, removedID[i]);

	return trimmedReads;
}

/*!
 * \brief Trims bad quality bases from two sets of sequences. This is done in a way
 * that ensures that the pair structure is conserved.
 * \param idSet1 StringSet of FastA-IDs for the first set of sequences.
 * \param seqSet1 StringSet containing the forward reads of the paired sequences.
 * \param idSet2 StringSet of FastA-IDs for the second set of sequences.
 * \param seqSet2 StringSet containing the backward reads of the paired sequences.
 * \param cutoff The minimum quality required from a base.
 * \param min_length The minimum length required after trimming. Shorter sequences
 *  will be deleted. (See remark).
 * \param spec The trimming algorithm used to trim the sequences.
 * \return A pair of unsigned ints containing the number of reads which had bases removed from.
 * \remark If one read is marked for removal, while its sibling is not, the read will
 * be replaced by a single "N". If both reads are marked for removal, they are removed.
 * This way the relation between paired reads stays intact.
 */
template <typename TId, typename TSeq, typename TSpec>
seqan::Pair<unsigned, unsigned> trimPairBatch(seqan::StringSet<TId> & idSet1, seqan::StringSet<TSeq> & seqSet1,
						  seqan::StringSet<TId> & idSet2, seqan::StringSet<TSeq> & seqSet2,
						  unsigned const cutoff, unsigned const min_length, TSpec const & spec)
{
	typedef typename seqan::Id<seqan::StringSet<TSeq> >::Type TSetID;
	seqan::String<TSetID> removedID1, removedID2;

	unsigned trimmedReads1 = _trimReads(removedID1, seqSet1, cutoff, min_length, spec);
	unsigned trimmedReads2 = _trimReads(removedID2, seqSet2, cutoff, min_length, spec);

	// Remove deleted sequences from both sets. If only one is removed, replace with N.
	// If both are removed, delete them both from the sets.
	int i = length(removedID1) - 1;
	int j = length(removedID2) - 1;

	// If one ID is larger than the largest of the paired set, it would
	// leave an orphaned paired read. If they both point to the same
	// position in their respective sets, we can remove both reads.
	while (i >= 0 && j >= 0)
	{
		if (removedID1[i] > removedID2[j])
			seqan::moveValue(seqSet1, removedID1[i--], TSeq("N"));
		else if (removedID2[j] > removedID1[i])
			seqan::moveValue(seqSet2, removedID2[j--], TSeq("N"));
		else
		{
			removeValues(idSet1, seqSet1, removedID1[i--]);
			removeValues(idSet2, seqSet2, removedID2[j--]);
		}
	}

	while (i >= 0) seqan::moveValue(seqSet1, removedID1[i--], TSeq("N"));
	while (j >= 0) seqan::moveValue(seqSet2, removedID2[j--], TSeq("N"));

	return seqan::Pair<unsigned, unsigned>(trimmedReads1, trimmedReads2);
}

#endif  // #ifndef SANDBOX_GROUP3_APPS_SEQDPT_READTRIMMING_H_
