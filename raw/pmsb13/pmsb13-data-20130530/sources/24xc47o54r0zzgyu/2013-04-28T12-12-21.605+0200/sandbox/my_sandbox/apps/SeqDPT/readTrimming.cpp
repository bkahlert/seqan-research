#include <iostream>
#include <seqan/sequence.h>  // CharString, ...
#include <seqan/seq_io.h>
#include <seqan/file.h>      // to stream a CharString into cout
#include "readTrimming.h"

// Trimming methods
// ----------------------------------------------------------------------------
// Simple trimming mechanism: Cut off bases as long as they are below our quality.
template <typename TSeq, typename TQual>
int _trimRead(TSeq& seq, TQual& qual, int const cutoff, int const min_length, Standard const &)
{
	int cut_pos = length(seq);
	int ret;
	for (int i=length(seq)-1; i >= min_length; i--)
		if ((seqan::ordValue(qual[i]) - 33) < cutoff)
			cut_pos = i;

	return cut_pos;
}

template <typename TSeq>
int _trimRead(TSeq& seq, int const cutoff, int const min_length, Standard const &)
{
	int cut_pos = length(seq);
	int ret;
	for (int i=length(seq)-1; i >= min_length; i--)
		if (seqan::getQualityValue(seq[i]) < cutoff)
			cut_pos = i;

	return cut_pos;
}

// Trimming mechanism used in BWA. Trim to argmax_x sum_{i=x+1}^l {cutoff - q_i}
template <typename TSeq, typename TQual>
int _trimRead(TSeq& seq, TQual& qual, int const cutoff, int const min_length, BWA const &)
{
	int max_arg = length(seq), sum = 0, max = 0;
	for (int i=length(seq)-1; i >= min_length; i--)
	{
		sum += cutoff - (seqan::ordValue(qual[i]) - 33);
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
int _trimRead(TSeq& seq, int const cutoff, int const min_length, BWA const &)
{
	int max_arg = length(seq), sum = 0, max = 0;
	for (int i=length(seq)-1; i >= min_length; i--)
	{
		sum += cutoff - seqan::getQualityValue(seq[i]);
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

// Trim by shifting a window over the text and calculating mean/median in it.
template <typename TSeq, typename TQual>
int _trimRead(TSeq& seq, TQual& qual, int const cutQual, int const min_length, Median const & spec)
{
	int window = spec.window;
	int bad_pos = length(seq);
	int cutoff = (cutQual+33)*window;

	// Shift window over read and keep mean quality.
	for (int i=min_length; i < length(seq); i++)
	{
		int avg = 0;
		// Take care only not to go over the end of the sequence. Shorten window near the end.
		int end = i + window < length(seq) ? window : length(seq) - i;

		// Calculate mean quality in the window.
		for (int j=0; j < end; j++)
			avg += seqan::ordValue(qual[i+j]);

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
			avg += seqan::getQualityValue(seq[i+j]);

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

// Trimming functions. With or without dedicated quality string.
// Return the number of bases trimmed off the end.
template <typename TSeq, typename TQual, typename TSpec>
int trimRead(TSeq& seq, TQual& qual, int const cutoff, int const min_length, TSpec const & spec)
{
	int ret, cut_pos;
	cut_pos = _trimRead(seq, qual, cutoff, min_length, spec);

	ret = length(seq) - cut_pos;
	seqan::erase(seq , cut_pos, length(seq));
	seqan::erase(qual, cut_pos, length(qual));

	return ret;
}

template <typename TSeq, typename TSpec>
int trimRead(TSeq& seq, int const cutoff, int const min_length, TSpec const & spec)
{
	int ret, cut_pos;
	cut_pos = _trimRead(seq, cutoff, min_length, spec);

	ret = length(seq) - cut_pos;
	seqan::erase(seq, cut_pos, length(seq));

	return ret;
}

int trimRead(seqan::String<seqan::Dna5Q>& seq, int const cutoff, int const min_length, Median const & spec)
{
	int ret, cut_pos;
	cut_pos = _trimRead(seq, cutoff, min_length, spec);

	ret = length(seq) - cut_pos;
	seqan::erase(seq, cut_pos, length(seq));

	return ret;
}

// trimBatch applies the trimRead function to each Sequence in a StringSet
// returns the number of trimmed sequences.
template <typename TSeq, typename TQual, typename TSpec>
int trimBatch(seqan::StringSet<TSeq> seqSet, seqan::StringSet<TQual> qualSet, int const cutoff, int const min_length, TSpec const & spec)
{
	int trimmedReads = 0;
	for (int i=0; i < length(seqSet); i++)
	{
		if (trimRead(value(seqSet, i), value(qualSet, i), cutoff, min_length, spec) > 0)
			trimmedReads++;
	}
	return trimmedReads;
}

template <typename TSeq, typename TSpec>
int trimBatch(seqan::StringSet<TSeq> seqSet, int const cutoff, int const min_length, TSpec const & spec)
{
	int trimmedReads = 0;
	for (int i=0; i < length(seqSet); i++)
	{
		if (trimRead(value(seqSet, i), cutoff, min_length, spec) > 0)
			trimmedReads++;
	}
	return trimmedReads;
}
