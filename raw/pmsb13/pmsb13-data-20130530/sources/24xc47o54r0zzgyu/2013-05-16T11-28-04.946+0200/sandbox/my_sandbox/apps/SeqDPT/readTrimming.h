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
	int window;
	Mean(int w) : window(w) {}
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
int getQuality(seqan::String<seqan::Dna5Q> seq, int i)
{
	return seqan::getQualityValue(seq[i]);
}

int getQuality(Dna5QAdapter seq, int i)
{
	return seqan::ordValue(seq.qual[i]) - 33;
}

// Length specialization for our Dna5Q wrapper struct.
int length(Dna5QAdapter seq)
{
	return length(seq.seq);
}

// Specialize seqan's erase function for our quality wrapper.
namespace seqan{
	void erase(Dna5QAdapter seq, int begin, int end)
	{
		erase(seq.seq, begin, end);
		erase(seq.qual, begin, end);
	}
}

// Trimming methods
// ----------------------------------------------------------------------------
// Simple trimming mechanism: Cut off bases as long as they are below our quality.
template <typename TSeq>
int _trimRead(TSeq& seq, int const cutoff, Tail const &)
{
	int cut_pos = length(seq);
	int ret;
	for (int i=length(seq)-1; i >= 0; i--)
		if (getQuality(seq, i) < cutoff)
			cut_pos = i;

	return cut_pos;
}

// Trimming mechanism used in BWA. Trim to argmax_x sum_{i=x+1}^l {cutoff - q_i}
template <typename TSeq>
int _trimRead(TSeq& seq, int const cutoff, BWA const &)
{
	int max_arg = length(seq), sum = 0, max = 0;
	for (int i=length(seq)-1; i >= 0; i--)
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
int _trimRead(TSeq& seq, int const cutQual, Mean const & spec)
{
	int window = spec.window;
	int avg = 0, i = 0;
	// Work with absolute cutoff in window to avoid divisions.
	int cutoff = cutQual*window;

	// Calculate average quality of initial window.
	for (i=0; i < window; i++)
		avg += getQuality(seq, i);

	// Shift window over read and keep mean quality, update in constant time.
	for (i=0; i < length(seq) && avg >= cutoff; i++)
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
int trimRead(TSeq& seq, TQual& qual, int const cutoff, TSpec const & spec)
{
	Dna5QAdapter a = Dna5QAdapter(seq, qual);
	return trimRead(a, cutoff, spec);
}

template <typename TSeq, typename TSpec>
int trimRead(TSeq& seq, int const cutoff, TSpec const & spec)
{
	int ret, cut_pos;
	cut_pos = _trimRead(seq, cutoff, spec);

	ret = length(seq) - cut_pos;
	seqan::erase(seq, cut_pos, length(seq));

	return ret;
}

// trimBatch applies the trimRead function to each Sequence in a StringSet
// returns the number of trimmed sequences.
template <typename TId, typename TSeq, typename TSpec>
int trimBatch(seqan::StringSet<TId> & idSet, seqan::StringSet<TSeq> & seqSet, int const cutoff, int const min_length, TSpec const & spec)
{
	int trimmedReads = 0;

	// Only save the StringSet-Ids of the strings that are going to be
	// removed. This avoids index problems when removing from the StringSet
	// in a running loop and enables easy parallelization with OpenMP.
	typedef typename seqan::Id<seqan::StringSet<TId> >::Type remID;
	typedef typename seqan::Id<seqan::StringSet<TSeq> >::Type remSeq;
	seqan::String<remID> removedID;
	seqan::String<remSeq> removedSeq;

	#pragma omp parallel for
	for (int i=0; i < length(seqSet); i++)
	{
		TSeq& read = value(seqSet, i);
		int trimmed = trimRead(read, cutoff, spec);
		if (trimmed > 0)
			#pragma omp atomic
			++trimmedReads;
		if (length(read) < min_length){
			#pragma omp critical
			{
				seqan::append(removedID, seqan::positionToId(idSet, i));
				seqan::append(removedSeq, seqan::positionToId(seqSet, i));
				//seqan::removeValueById(seqSet, seqan::positionToId(seqSet, i));
				//seqan::removeValueById(idSet, seqan::positionToId(idSet, i));
				//--i; // Decrease index to compensate for the shift in positions caused by the deletion.
			}
		}
	}

	for (int i=0; i < length(removedID); i++)
	{
		seqan::removeValueById(idSet, removedID[i]);
		seqan::removeValueById(seqSet, removedSeq[i]);
	}


	return trimmedReads;
}

#endif  // #ifndef SANDBOX_GROUP3_APPS_SEQDPT_READTRIMMING_H_
