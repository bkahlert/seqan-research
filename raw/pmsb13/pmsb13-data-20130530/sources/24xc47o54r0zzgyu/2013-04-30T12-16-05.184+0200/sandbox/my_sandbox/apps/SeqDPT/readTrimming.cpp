/*#include <iostream>
#include <seqan/sequence.h>  // CharString, ...
#include <seqan/seq_io.h>
#include <seqan/file.h>      // to stream a CharString into cout
#include "readTrimming.h"

struct Dna5QAdapter{
	seqan::String<seqan::Dna5> seq;
	seqan::CharString qual;
	Dna5QAdapter(seqan::String<seqan::Dna5> s, seqan::CharString q)
	{
		seq = s;
		qual = q;
	}
};

int getQuality(seqan::String<seqan::Dna5Q> seq, int i)
{
	return seqan::getQualityValue(seq[i]);
}

int getQuality(Dna5QAdapter seq, int i)
{
	return seqan::ordValue(seq.qual[i]) - 33;
}

// Trimming methods
// ----------------------------------------------------------------------------
// Simple trimming mechanism: Cut off bases as long as they are below our quality.
template <typename TSeq>
int _trimRead(TSeq& seq, int const cutoff, int const min_length, Standard const &)
{
	int cut_pos = length(seq);
	int ret;
	for (int i=length(seq)-1; i >= min_length; i--)
		if (getQuality(seq, i) < cutoff)
			cut_pos = i;

	return cut_pos;
}

// Trimming mechanism used in BWA. Trim to argmax_x sum_{i=x+1}^l {cutoff - q_i}
template <typename TSeq>
int _trimRead(TSeq& seq, int const cutoff, int const min_length, BWA const &)
{
	int max_arg = length(seq), sum = 0, max = 0;
	for (int i=length(seq)-1; i >= min_length; i--)
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

template <typename TSeq>
int _trimRead(TSeq& seq, int const cutQual, int const min_length, Median const & spec)
{
	int window = spec.window;
	int bad_pos = length(seq);
	int cutoff = cutQual*window;

	// Shift window over read and keep mean quality.
	for (int i=min_length; i < length(seq); i++)
	{
		int avg = 0;
		// Take care only not to go over the end of the sequence. Shorten window near the end.
		int end = i + window < length(seq) ? window : length(seq) - i;

		// Calculate mean quality in the window.
		for (int j=0; j < end; j++)
			avg += getQuality(seq, i+j);

		// If this window was bad, we stop here and cut the tail.
		if (avg < cutoff)
		{
			bad_pos = i;
			break;
		}
	}

	// We will cut if we found a bad window <=> if bad_pos was changed.
	return bad_pos;
}
*/
