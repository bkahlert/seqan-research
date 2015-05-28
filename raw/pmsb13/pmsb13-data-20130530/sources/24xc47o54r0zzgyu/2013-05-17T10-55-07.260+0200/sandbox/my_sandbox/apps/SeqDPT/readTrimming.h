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
// Author: Your Name <your.email@example.net>
// ==========================================================================

#ifndef SANDBOX_GROUP3_APPS_SEQDPT_READTRIMMING_H_
#define SANDBOX_GROUP3_APPS_SEQDPT_READTRIMMING_H_

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// Tags for choosing the appropriate trimming method.
struct TrimmingAlgorithm {};

struct Tail : TrimmingAlgorithm {};
struct BWA : TrimmingAlgorithm {};
struct Mean : TrimmingAlgorithm {
	unsigned window;
	Mean(unsigned w) : window(w) {}
};

struct TrimConfig{
	unsigned cutoff, min_length;
	TrimConfig(unsigned c, unsigned m) : cutoff(c), min_length(m) {}
};

// Wraps sequence and quality, making them usable for our trimRead function.
struct Dna5QAdapter{
	seqan::String<seqan::Dna5>& seq;
	seqan::CharString& qual;
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

// Get the phred quality of a sequence position.
unsigned getQuality(const seqan::String<seqan::Dna5Q>& seq, unsigned i)
{
	return seqan::getQualityValue(seq[i]);
}

unsigned getQuality(const Dna5QAdapter& seq, unsigned i)
{
	return seqan::ordValue(seq.qual[i]) - 33;
}


// Length specialization for our Dna5Q wrapper struct.
unsigned length(const Dna5QAdapter & seq)
{
	return length(seq.seq);
}

// Specialize seqan's erase function for our quality wrapper.
void erase(const Dna5QAdapter & seq, unsigned begin, unsigned end)
{
	seqan::erase(seq.seq , begin, end);
	seqan::erase(seq.qual, begin, end);
}

// Trimming methods
// ----------------------------------------------------------------------------
// Simple trimming mechanism: Cut off bases as long as they are below our quality.
template <typename TSeq>
unsigned _trimRead(const TSeq& seq, unsigned const cutoff, Tail const &)
{
	unsigned cut_pos = length(seq);
	unsigned ret;
	for (int i=length(seq)-1; i >= 0; --i)
		if (getQuality(seq, i) < cutoff)
			cut_pos = i;

	return cut_pos;
}

// Trimming mechanism used in BWA. Trim to argmax_x sum_{i=x+1}^l {cutoff - q_i}
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

// Trimming functions. With or without dedicated quality string.
// Return the number of bases trimmed off the end.
template <typename TSeq, typename TQual, typename TSpec>
unsigned trimRead(TSeq& seq, TQual& qual, unsigned const cutoff, TSpec const & spec)
{
	Dna5QAdapter a = Dna5QAdapter(seq, qual);
	return trimRead(a, cutoff, spec);
}

template <typename TSeq, typename TSpec>
unsigned trimRead(TSeq& seq, unsigned const cutoff, TSpec const & spec)
{
	unsigned ret, cut_pos;
	cut_pos = _trimRead(seq, cutoff, spec);

	ret = length(seq) - cut_pos;
	erase(seq, cut_pos, length(seq));

	return ret;
}

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

template <typename TSet1, typename TSet2, typename TId>
void removeValues(TSet1 & set1, TSet2 & set2, TId id)
{
	seqan::removeValueById( set1, id);
	seqan::removeValueById( set2, id);
}

/**
 * trimBatch applies the trimRead function to each Sequence in a StringSet
 * returns the number of trimmed sequences.
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

/**
 * trimPairBatch does the same as trimBatch, but keeps track of paired reads
 * and in case of removed reads either replaces missing reads with an N or
 * an the case of a fully removed pair, deletes both paired reads.
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
